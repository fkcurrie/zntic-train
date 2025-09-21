import os
import pandas as pd
from google.cloud import storage
import io
import joblib
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
import numpy as np
import json
import argparse
import tensorflow as tf

def create_neural_network(input_shape, hyperparams):
    """Builds and compiles a Keras Sequential model."""
    layers = hyperparams.get("layers", [128, 64])
    activation = hyperparams.get("activation", "relu")
    optimizer = hyperparams.get("optimizer", "adam")
    
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Input(shape=(input_shape,)))
    
    for units in layers:
        model.add(tf.keras.layers.Dense(units, activation=activation))
        
    model.add(tf.keras.layers.Dense(1, activation='sigmoid'))
    
    model.compile(optimizer=optimizer,
                  loss='binary_crossentropy',
                  metrics=['accuracy'])
    return model

def main(args):
    """
    Main function to train the model.
    """
    storage_client = storage.Client()
    bucket = storage_client.get_bucket("zntic-data")

    print("Downloading features from GCS...")
    blob = bucket.blob("features.csv")
    features_data = blob.download_as_string()

    print("Loading features into DataFrame...")
    df = pd.read_csv(io.StringIO(features_data.decode("utf-8")))

    print("Generating placeholder labels...")
    np.random.seed(42)
    df['zoonotic'] = np.random.randint(0, 2, df.shape[0])

    print("Splitting data into training and testing sets...")
    X = df.drop('zoonotic', axis=1)
    y = df['zoonotic']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    print(f"Saving feature columns to GCS for model '{args.model_name}'...")
    feature_columns = X_train.columns.tolist()
    columns_path = f"models/{args.model_name}/feature_columns.json"
    columns_blob = bucket.blob(columns_path)
    columns_blob.upload_from_string(json.dumps(feature_columns))
    
    print(f"Training {args.model_type} model...")
    hyperparams = json.loads(args.hyperparams)
    
    if args.model_type == 'NeuralNetwork':
        model = create_neural_network(X_train.shape[1], hyperparams)
        epochs = hyperparams.get("epochs", 10)
        batch_size = hyperparams.get("batch_size", 32)
        model.fit(X_train, y_train, epochs=epochs, batch_size=batch_size, validation_split=0.2)
        # For evaluation, we need to get the raw predictions
        y_pred_proba = model.predict(X_test)
        y_pred = (y_pred_proba > 0.5).astype(int)
    else:
        if args.model_type == 'LogisticRegression':
            model = LogisticRegression(**hyperparams)
        elif args.model_type == 'RandomForestClassifier':
            model = RandomForestClassifier(**hyperparams)
        elif args.model_type == 'SVC':
            hyperparams['probability'] = True
            model = SVC(**hyperparams)
        elif args.model_type == 'GradientBoostingClassifier':
            model = GradientBoostingClassifier(**hyperparams)
        else:
            raise ValueError(f"Unsupported model type: {args.model_type}")
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

    acc = accuracy_score(y_test, y_pred)
    print(f"Model Accuracy: {acc:.2f}")

    print("Saving model to local file...")
    if args.model_type == 'NeuralNetwork':
        model_filename = 'model.keras'
        model.save(model_filename)
    else:
        model_filename = 'model.joblib'
        joblib.dump(model, model_filename)

    print(f"Uploading model to GCS at 'models/{args.model_name}/'...")
    model_path = f"models/{args.model_name}/{model_filename}"
    model_blob = bucket.blob(model_path)
    model_blob.upload_from_filename(model_filename)

    print(f"Training complete. Model '{args.model_name}' and feature columns uploaded to GCS.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a zoonotic potential prediction model.")
    parser.add_argument(
        '--model-type', 
        type=str, 
        required=True, 
        choices=['LogisticRegression', 'RandomForestClassifier', 'SVC', 'GradientBoostingClassifier', 'NeuralNetwork'],
        help='The type of model to train.'
    )
    parser.add_argument(
        '--model-name', 
        type=str, 
        required=True, 
        help='The name to save the model under in GCS.'
    )
    parser.add_argument(
        '--hyperparams', 
        type=str, 
        default='{}', 
        help='A JSON string of hyperparameters for the model.'
    )
    
    parsed_args = parser.parse_args()
    main(parsed_args)
