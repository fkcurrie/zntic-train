class RetrainingDashboard {
    constructor() {
        this.apiUrl = 'http://34.30.222.113'; // Assuming the same API URL
        this.init();
    }

    init() {
        this.setupEventListeners();
        this.updateHyperparameterPlaceholder();
    }

    setupEventListeners() {
        document.getElementById('model-type').addEventListener('change', () => {
            this.updateHyperparameterPlaceholder();
        });

        document.getElementById('start-training-btn').addEventListener('click', () => {
            this.startTraining();
        });
    }

    updateHyperparameterPlaceholder() {
        const modelType = document.getElementById('model-type').value;
        const textarea = document.getElementById('hyperparams');
        let placeholder = {};

        if (modelType === 'RandomForestClassifier') {
            placeholder = {
                "n_estimators": 100,
                "max_depth": 10,
                "random_state": 42
            };
        } else if (modelType === 'LogisticRegression') {
            placeholder = {
                "C": 1.0,
                "max_iter": 1000,
                "solver": "lbfgs"
            };
        } else if (modelType === 'SVC') {
            placeholder = {
                "C": 1.0,
                "kernel": "rbf",
                "gamma": "scale",
                "random_state": 42
            };
        } else if (modelType === 'GradientBoostingClassifier') {
            placeholder = {
                "n_estimators": 100,
                "learning_rate": 0.1,
                "max_depth": 3,
                "random_state": 42
            };
        } else if (modelType === 'NeuralNetwork') {
            placeholder = {
                "layers": [128, 64],
                "activation": "relu",
                "optimizer": "adam",
                "epochs": 50,
                "batch_size": 32
            };
        }


        textarea.placeholder = JSON.stringify(placeholder, null, 4);
        // Clear the user's text to show the new placeholder
        textarea.value = '';
    }

    async startTraining() {
        this.showLoading();
        this.hideStatus();
        this.hideError();

        const modelName = document.getElementById('model-name').value.trim();
        const modelType = document.getElementById('model-type').value;
        const hyperparamsRaw = document.getElementById('hyperparams').value.trim();

        if (!modelName) {
            this.showError('Model Name is required.');
            this.hideLoading();
            return;
        }

        let hyperparams;
        try {
            // If the user leaves the textarea blank, use the placeholder.
            if (hyperparamsRaw) {
                hyperparams = JSON.parse(hyperparamsRaw);
            } else {
                hyperparams = JSON.parse(document.getElementById('hyperparams').placeholder);
            }
        } catch (e) {
            this.showError('Invalid JSON in Hyperparameters field.');
            this.hideLoading();
            return;
        }

        const payload = {
            model_name: modelName,
            model_type: modelType,
            hyperparams: hyperparams
        };

        try {
            const response = await fetch(`${this.apiUrl}/retrain`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(payload),
            });

            const result = await response.json();

            if (!response.ok) {
                throw new Error(result.error || `HTTP error! status: ${response.status}`);
            }

            this.showStatus(`Success! Training job submitted as '${result.job_name}'. You can now go back to the predictor page. The new model will appear in the dropdown once training is complete.`);

        } catch (error) {
            console.error('Training submission failed:', error);
            this.showError(`Submission failed: ${error.message}`);
        } finally {
            this.hideLoading();
        }
    }

    showLoading() {
        document.getElementById('loading').classList.remove('hidden');
        document.getElementById('start-training-btn').disabled = true;
    }

    hideLoading() {
        document.getElementById('loading').classList.add('hidden');
        document.getElementById('start-training-btn').disabled = false;
    }

    showStatus(message) {
        document.getElementById('status-message').textContent = message;
        document.getElementById('status').classList.remove('hidden');
    }

    hideStatus() {
        document.getElementById('status').classList.add('hidden');
    }

    showError(message) {
        document.getElementById('error-message').textContent = message;
        document.getElementById('error').classList.remove('hidden');
    }

    hideError() {
        document.getElementById('error').classList.add('hidden');
    }
}

document.addEventListener('DOMContentLoaded', () => {
    new RetrainingDashboard();
});
