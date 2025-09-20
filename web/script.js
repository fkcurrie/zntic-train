// Modern JavaScript for Zoonotic Predictor
class ZoonoticPredictor {
    constructor() {
        this.apiUrl = 'http://34.30.222.113';
        this.currentTab = 'manual';
        this.init();
    }

    init() {
        this.setupEventListeners();
        this.checkApiStatus();
        this.setupTabSystem();
        this.fetchBuildDates();
        this.fetchModels();
    }

    setupEventListeners() {
        // Tab switching
        document.querySelectorAll('.tab-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                this.switchTab(e.target.dataset.tab);
            });
        });

        // Sample data buttons
        document.querySelectorAll('.sample-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const button = e.target.closest('.sample-btn');
                if (button && button.dataset.sample) {
                    this.loadSampleData(button.dataset.sample);
                    setTimeout(() => this.analyzeData(), 100);
                }
            });
        });

        // Control buttons
        document.getElementById('analyze-btn').addEventListener('click', () => {
            this.analyzeData();
        });

        document.getElementById('clear-btn').addEventListener('click', () => {
            this.clearData();
        });
        
        document.getElementById('refresh-models-btn').addEventListener('click', () => {
            this.fetchModels();
        });

        // Auto-resize textareas
        document.querySelectorAll('textarea').forEach(textarea => {
            textarea.addEventListener('input', () => {
                this.autoResizeTextarea(textarea);
            });
        });
    }

    setupTabSystem() {
        this.switchTab('manual');
    }

    switchTab(tabName) {
        this.currentTab = tabName;
        document.querySelectorAll('.tab-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.tab === tabName);
        });
        document.querySelectorAll('.tab-content').forEach(content => {
            content.classList.toggle('active', content.id === `${tabName}-tab`);
        });
        this.hideResults();
        this.hideError();
    }

    async fetchModels() {
        const selector = document.getElementById('model-selector');
        selector.innerHTML = '<option>Loading models...</option>';
        selector.disabled = true;

        try {
            const response = await fetch(`${this.apiUrl}/models`);
            if (!response.ok) throw new Error('Failed to fetch models');
            
            const models = await response.json();
            selector.innerHTML = ''; // Clear loading text
            
            if (models.length === 0) {
                selector.innerHTML = '<option>No models found</option>';
            } else {
                models.forEach(modelName => {
                    const option = document.createElement('option');
                    option.value = modelName;
                    option.textContent = modelName;
                    selector.appendChild(option);
                });
            }
        } catch (error) {
            console.error('Failed to fetch models:', error);
            selector.innerHTML = '<option>Error loading models</option>';
        } finally {
            selector.disabled = false;
        }
    }

    async checkApiStatus() {
        const statusElement = document.getElementById('api-status');
        const statusIndicator = statusElement.querySelector('.status-indicator');
        const statusText = statusElement.querySelector('span');

        try {
            const controller = new AbortController();
            const timeoutId = setTimeout(() => controller.abort(), 5000);
            const response = await fetch(`${this.apiUrl}/health`, { signal: controller.signal });
            clearTimeout(timeoutId);

            if (response.ok) {
                statusIndicator.classList.add('online');
                statusIndicator.classList.remove('offline');
                statusText.textContent = 'API Online';
            } else {
                throw new Error('API not responding correctly');
            }
        } catch (error) {
            statusIndicator.classList.add('offline');
            statusIndicator.classList.remove('online');
            statusText.textContent = 'API Offline';
            console.error('API status check failed:', error);
        }
    }

    loadSampleData(sampleType) {
        let features;
        switch (sampleType) {
            case 'high-risk':
                features = this.generateHighRiskSample();
                break;
            case 'low-risk':
                features = this.generateLowRiskSample();
                break;
            case 'random':
                features = this.generateRandomSample();
                break;
            default:
                return;
        }
        this.switchTab('manual');
        const featuresInput = document.getElementById('features-input');
        featuresInput.value = features.join(', ');
        this.autoResizeTextarea(featuresInput);
    }

    generateHighRiskSample() {
        return Array.from({ length: 64 }, () => parseFloat((Math.random() * 0.4 + 0.6).toFixed(4)));
    }

    generateLowRiskSample() {
        return Array.from({ length: 64 }, () => parseFloat((Math.random() * 0.05).toFixed(4)));
    }

    generateRandomSample() {
        return Array.from({ length: 64 }, () => parseFloat(Math.random().toFixed(4)));
    }

    async analyzeData() {
        this.showLoading();
        this.hideResults();
        this.hideError();

        const modelName = document.getElementById('model-selector').value;
        if (!modelName || modelName === 'No models found' || modelName === 'Error loading models') {
            this.showError('Please select a valid model.');
            this.hideLoading();
            return;
        }

        const features = this.extractFeatures();
        if (!features) {
            this.showError('Please provide valid genomic features.');
            this.hideLoading();
            return;
        }

        try {
            const controller = new AbortController();
            const timeoutId = setTimeout(() => controller.abort(), 30000);
            const response = await fetch(`${this.apiUrl}/predict`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ model_name: modelName, features: features }),
                signal: controller.signal
            });
            clearTimeout(timeoutId);

            const result = await response.json();
            if (!response.ok) {
                throw new Error(result.error || `HTTP error! status: ${response.status}`);
            }
            this.displayResults(result);
        } catch (error) {
            console.error('Analysis failed:', error);
            this.showError(`Analysis failed: ${error.message}`);
        } finally {
            this.hideLoading();
        }
    }

    extractFeatures() {
        const tab = this.currentTab;
        if (tab === 'manual' || tab === 'sample') {
            return this.extractManualFeatures();
        } else if (tab === 'sequence') {
            return this.extractSequenceFeatures();
        }
        return null;
    }

    extractManualFeatures() {
        const input = document.getElementById('features-input').value.trim();
        if (!input) return null;
        try {
            const features = input.split(',').map(val => parseFloat(val.trim())).filter(val => !isNaN(val));
            if (features.length !== 64) throw new Error(`Expected 64 features, got ${features.length}`);
            return features;
        } catch (error) {
            this.showError(`Invalid feature format: ${error.message}`);
            return null;
        }
    }

    extractSequenceFeatures() {
        const sequence = document.getElementById('sequence-input').value.trim();
        if (!sequence) return null;
        this.showError('Sequence analysis feature is not yet implemented. Please use manual input or sample data.');
        return null;
    }

    displayResults(result) {
        const predictionBadge = document.getElementById('prediction-badge');
        const predictionIcon = document.getElementById('prediction-icon');
        const predictionText = document.getElementById('prediction-text');

        predictionBadge.className = `prediction-badge ${result.prediction}`;
        if (result.prediction === 'zoonotic') {
            predictionIcon.className = 'fas fa-exclamation-triangle';
            predictionText.textContent = 'High Zoonotic Potential';
        } else {
            predictionIcon.className = 'fas fa-check-circle';
            predictionText.textContent = 'Low Zoonotic Potential';
        }

        const nonZoonoticPercent = (result.confidence['non-zoonotic'] * 100).toFixed(1);
        const zoonoticPercent = (result.confidence.zoonotic * 100).toFixed(1);

        document.getElementById('non-zoonotic-bar').style.width = `${nonZoonoticPercent}%`;
        document.getElementById('zoonotic-bar').style.width = `${zoonoticPercent}%`;
        document.getElementById('non-zoonotic-percent').textContent = `${nonZoonoticPercent}%`;
        document.getElementById('zoonotic-percent').textContent = `${zoonoticPercent}%`;
        document.getElementById('interpretation-text').textContent = this.generateInterpretation(result);
        this.showResults();
    }

    generateInterpretation(result) {
        const zoonoticProb = result.confidence.zoonotic;
        if (result.prediction === 'zoonotic') {
            if (zoonoticProb > 0.8) return `The model predicts a HIGH zoonotic potential (${(zoonoticProb * 100).toFixed(1)}% confidence). This viral genome shows strong patterns associated with viruses capable of human transmission.`;
            if (zoonoticProb > 0.6) return `The model predicts MODERATE zoonotic potential (${(zoonoticProb * 100).toFixed(1)}% confidence). While concerning, this requires further investigation and monitoring.`;
            return `The model suggests possible zoonotic potential (${(zoonoticProb * 100).toFixed(1)}% confidence), but with lower certainty. Continue monitoring.`;
        } else {
            if (zoonoticProb < 0.2) return `The model indicates LOW zoonotic potential (${((1 - zoonoticProb) * 100).toFixed(1)}% non-zoonotic confidence). This viral genome shows patterns typical of viruses that do not readily transmit to humans.`;
            if (zoonoticProb < 0.4) return `The model suggests relatively low zoonotic risk (${((1 - zoonoticProb) * 100).toFixed(1)}% non-zoonotic confidence). Standard monitoring protocols should be sufficient.`;
            return `While classified as non-zoonotic, the confidence is moderate (${((1 - zoonoticProb) * 100).toFixed(1)}%). Consider additional analysis to confirm the assessment.`;
        }
    }

    clearData() {
        document.getElementById('features-input').value = '';
        document.getElementById('sequence-input').value = '';
        document.querySelectorAll('textarea').forEach(textarea => {
            textarea.style.height = 'auto';
        });
        this.hideResults();
        this.hideError();
    }

    autoResizeTextarea(textarea) {
        textarea.style.height = 'auto';
        textarea.style.height = `${textarea.scrollHeight}px`;
    }

    showLoading() {
        document.getElementById('loading').classList.remove('hidden');
        document.getElementById('analyze-btn').disabled = true;
    }

    hideLoading() {
        document.getElementById('loading').classList.add('hidden');
        document.getElementById('analyze-btn').disabled = false;
    }

    showResults() {
        const resultsEl = document.getElementById('results');
        resultsEl.classList.remove('hidden');
        resultsEl.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }

    hideResults() {
        document.getElementById('results').classList.add('hidden');
    }

    showError(message) {
        const errorEl = document.getElementById('error');
        document.getElementById('error-message').textContent = message;
        errorEl.classList.remove('hidden');
        errorEl.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }

    hideError() {
        document.getElementById('error').classList.add('hidden');
    }

    async fetchBuildDates() {
        try {
            const response = await fetch('/build-info.json');
            const info = await response.json();
            document.getElementById('web-build-date').textContent = info.build_date || 'Static';
        } catch (error) {
            console.error('Failed to fetch web build date:', error);
            document.getElementById('web-build-date').textContent = 'Static';
        }

        try {
            const controller = new AbortController();
            const timeoutId = setTimeout(() => controller.abort(), 5000);
            const response = await fetch(`${this.apiUrl}/info`, { signal: controller.signal });
            clearTimeout(timeoutId);
            const info = await response.json();
            document.getElementById('model-build-date').textContent = info.model_training_date || 'Unknown';
        } catch (error) {
            console.error('Failed to fetch model build date:', error);
            document.getElementById('model-build-date').textContent = 'Unavailable';
        }
    }
}

document.addEventListener('DOMContentLoaded', () => {
    new ZoonoticPredictor();
});