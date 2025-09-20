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
                // Find the button element (in case user clicks on icon or text inside)
                const button = e.target.closest('.sample-btn');
                if (button && button.dataset.sample) {
                    this.loadSampleData(button.dataset.sample);
                    // Automatically run analysis after loading sample data
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

        // Auto-resize textareas
        document.querySelectorAll('textarea').forEach(textarea => {
            textarea.addEventListener('input', () => {
                this.autoResizeTextarea(textarea);
            });
        });
    }

    setupTabSystem() {
        // Ensure proper initial state
        this.switchTab('manual');
    }

    switchTab(tabName) {
        this.currentTab = tabName;

        // Update tab buttons
        document.querySelectorAll('.tab-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.tab === tabName);
        });

        // Update tab content
        document.querySelectorAll('.tab-content').forEach(content => {
            content.classList.toggle('active', content.id === `${tabName}-tab`);
        });

        // Clear any existing data when switching tabs
        this.hideResults();
        this.hideError();
    }

    async checkApiStatus() {
        const statusElement = document.getElementById('api-status');
        const statusIndicator = statusElement.querySelector('.status-indicator');
        const statusText = statusElement.querySelector('span');

        try {
            const controller = new AbortController();
            const timeoutId = setTimeout(() => controller.abort(), 5000);

            const response = await fetch(`${this.apiUrl}/health`, {
                method: 'GET',
                signal: controller.signal
            });

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
                // Generate sample data that would likely be classified as zoonotic
                features = this.generateHighRiskSample();
                break;
            case 'low-risk':
                // Generate sample data that would likely be classified as non-zoonotic
                features = this.generateLowRiskSample();
                break;
            case 'random':
                // Generate completely random sample data
                features = this.generateRandomSample();
                break;
            default:
                return;
        }

        // Switch to manual tab and populate the textarea
        this.switchTab('manual');
        document.getElementById('features-input').value = features.join(', ');

        // Auto-resize the textarea
        this.autoResizeTextarea(document.getElementById('features-input'));
    }

    generateHighRiskSample() {
        // Create sample data with higher values that might indicate zoonotic potential
        const features = [];
        for (let i = 0; i < 64; i++) {
            // Bias towards higher values for "high risk"
            const value = Math.random() * 0.4 + 0.6; // 0.6 to 1.0
            features.push(parseFloat(value.toFixed(4)));
        }
        return features;
    }

    generateLowRiskSample() {
        // Create sample data with lower values that might indicate non-zoonotic
        const features = [];
        for (let i = 0; i < 64; i++) {
            // Much lower values for "low risk" - very conservative range
            const value = Math.random() * 0.05; // 0.0 to 0.05 (extremely low range)
            features.push(parseFloat(value.toFixed(4)));
        }
        return features;
    }

    generateRandomSample() {
        // Create completely random sample data
        const features = [];
        for (let i = 0; i < 64; i++) {
            const value = Math.random(); // 0.0 to 1.0
            features.push(parseFloat(value.toFixed(4)));
        }
        return features;
    }

    async analyzeData() {
        try {
            this.showLoading();
            this.hideResults();
            this.hideError();

            const features = this.extractFeatures();
            if (!features) {
                this.showError('Please provide valid genomic features');
                return;
            }

            const controller = new AbortController();
            const timeoutId = setTimeout(() => controller.abort(), 30000); // 30 second timeout for ML predictions

            const response = await fetch(`${this.apiUrl}/predict`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({ features: features }),
                signal: controller.signal
            });

            clearTimeout(timeoutId);

            if (!response.ok) {
                const errorData = await response.json();
                throw new Error(errorData.error || `HTTP error! status: ${response.status}`);
            }

            const result = await response.json();
            this.displayResults(result);

        } catch (error) {
            console.error('Analysis failed:', error);
            this.showError(`Analysis failed: ${error.message}`);
        } finally {
            this.hideLoading();
        }
    }

    extractFeatures() {
        switch (this.currentTab) {
            case 'manual':
                return this.extractManualFeatures();
            case 'sequence':
                return this.extractSequenceFeatures();
            case 'sample':
                return this.extractManualFeatures(); // Sample data goes to manual input
            default:
                return null;
        }
    }

    extractManualFeatures() {
        const input = document.getElementById('features-input').value.trim();
        if (!input) {
            return null;
        }

        try {
            // Split by comma and clean up the values
            const features = input.split(',')
                .map(val => parseFloat(val.trim()))
                .filter(val => !isNaN(val));

            if (features.length !== 64) {
                throw new Error(`Expected 64 features, got ${features.length}`);
            }

            return features;
        } catch (error) {
            this.showError(`Invalid feature format: ${error.message}`);
            return null;
        }
    }

    extractSequenceFeatures() {
        const sequence = document.getElementById('sequence-input').value.trim();
        if (!sequence) {
            return null;
        }

        // This is a simplified version - in reality, you'd need sophisticated
        // bioinformatics tools to extract features from genomic sequences
        this.showError('Sequence analysis feature is not yet implemented. Please use manual input or sample data.');
        return null;
    }

    displayResults(result) {
        const resultsElement = document.getElementById('results');

        // Update prediction badge
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

        // Update confidence bars - handle different API response formats
        let nonZoonoticPercent, zoonoticPercent;

        if (result.confidence && typeof result.confidence.non_zoonotic !== 'undefined') {
            nonZoonoticPercent = (result.confidence.non_zoonotic * 100).toFixed(1);
            zoonoticPercent = (result.confidence.zoonotic * 100).toFixed(1);
        } else if (result.confidence && Array.isArray(result.confidence)) {
            // Handle array format [non_zoonotic, zoonotic]
            nonZoonoticPercent = (result.confidence[0] * 100).toFixed(1);
            zoonoticPercent = (result.confidence[1] * 100).toFixed(1);
        } else {
            // Fallback - calculate from prediction
            if (result.prediction === 'zoonotic') {
                zoonoticPercent = '75.0';
                nonZoonoticPercent = '25.0';
            } else {
                nonZoonoticPercent = '75.0';
                zoonoticPercent = '25.0';
            }
        }

        document.getElementById('non-zoonotic-bar').style.width = `${nonZoonoticPercent}%`;
        document.getElementById('zoonotic-bar').style.width = `${zoonoticPercent}%`;
        document.getElementById('non-zoonotic-percent').textContent = `${nonZoonoticPercent}%`;
        document.getElementById('zoonotic-percent').textContent = `${zoonoticPercent}%`;

        // Update interpretation
        const interpretationText = this.generateInterpretation(result);
        document.getElementById('interpretation-text').textContent = interpretationText;

        // Show results with animation
        this.showResults();
    }

    generateInterpretation(result) {
        const confidence = result.confidence;
        let zoonoticProb, nonZoonoticProb;

        // Handle different confidence formats
        if (confidence && typeof confidence.zoonotic !== 'undefined') {
            zoonoticProb = confidence.zoonotic;
            nonZoonoticProb = confidence.non_zoonotic;
        } else if (confidence && Array.isArray(confidence)) {
            nonZoonoticProb = confidence[0];
            zoonoticProb = confidence[1];
        } else {
            // Fallback values
            zoonoticProb = result.prediction === 'zoonotic' ? 0.75 : 0.25;
            nonZoonoticProb = result.prediction === 'zoonotic' ? 0.25 : 0.75;
        }

        if (result.prediction === 'zoonotic') {
            if (zoonoticProb > 0.8) {
                return `The model predicts a HIGH zoonotic potential (${(zoonoticProb * 100).toFixed(1)}% confidence). This viral genome shows strong patterns associated with viruses capable of human transmission. Immediate surveillance and containment measures may be warranted.`;
            } else if (zoonoticProb > 0.6) {
                return `The model predicts MODERATE zoonotic potential (${(zoonoticProb * 100).toFixed(1)}% confidence). While concerning, this requires further investigation and monitoring. Consider enhanced biosecurity measures.`;
            } else {
                return `The model suggests possible zoonotic potential (${(zoonoticProb * 100).toFixed(1)}% confidence), but with lower certainty. Continue monitoring and consider additional testing.`;
            }
        } else {
            if (zoonoticProb < 0.2) {
                return `The model indicates LOW zoonotic potential (${(nonZoonoticProb * 100).toFixed(1)}% non-zoonotic confidence). This viral genome shows patterns typical of viruses that do not readily transmit to humans.`;
            } else if (zoonoticProb < 0.4) {
                return `The model suggests relatively low zoonotic risk (${(nonZoonoticProb * 100).toFixed(1)}% non-zoonotic confidence). Standard monitoring protocols should be sufficient.`;
            } else {
                return `While classified as non-zoonotic, the confidence is moderate (${(nonZoonoticProb * 100).toFixed(1)}%). Consider additional analysis to confirm the assessment.`;
            }
        }
    }

    clearData() {
        // Clear all input fields
        document.getElementById('features-input').value = '';
        document.getElementById('sequence-input').value = '';

        // Reset textareas to normal size
        document.querySelectorAll('textarea').forEach(textarea => {
            textarea.style.height = 'auto';
        });

        // Hide results and errors
        this.hideResults();
        this.hideError();
    }

    autoResizeTextarea(textarea) {
        textarea.style.height = 'auto';
        textarea.style.height = textarea.scrollHeight + 'px';
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
        document.getElementById('results').classList.remove('hidden');
        // Smooth scroll to results
        document.getElementById('results').scrollIntoView({
            behavior: 'smooth',
            block: 'start'
        });
    }

    hideResults() {
        document.getElementById('results').classList.add('hidden');
    }

    showError(message) {
        const errorElement = document.getElementById('error');
        document.getElementById('error-message').textContent = message;
        errorElement.classList.remove('hidden');

        // Smooth scroll to error
        errorElement.scrollIntoView({
            behavior: 'smooth',
            block: 'start'
        });
    }

    hideError() {
        document.getElementById('error').classList.add('hidden');
    }

    async fetchBuildDates() {
        // Set web interface build date to current time (since this deployment is happening now)
        const webBuildDate = new Date().toISOString().slice(0, 19).replace('T', ' ') + ' UTC';
        document.getElementById('web-build-date').textContent = webBuildDate;

        // Try to fetch API build information
        try {
            const controller = new AbortController();
            const timeoutId = setTimeout(() => controller.abort(), 5000);

            const response = await fetch(`${this.apiUrl}/info`, {
                method: 'GET',
                signal: controller.signal
            });

            clearTimeout(timeoutId);

            if (response.ok) {
                const info = await response.json();

                // Display API build date if available
                if (info.build_date) {
                    document.getElementById('api-build-date').textContent = info.build_date;
                } else {
                    document.getElementById('api-build-date').textContent = 'Unknown';
                }

                // Display model training date if available
                if (info.model_training_date) {
                    document.getElementById('model-build-date').textContent = info.model_training_date;
                } else {
                    document.getElementById('model-build-date').textContent = 'Unknown';
                }
            } else {
                throw new Error('API info not available');
            }
        } catch (error) {
            console.error('Failed to fetch build dates:', error);
            document.getElementById('api-build-date').textContent = 'Unavailable';
            document.getElementById('model-build-date').textContent = 'Unavailable';
        }
    }
}

// Initialize the application when the DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new ZoonoticPredictor();
});

// Add some visual enhancements
document.addEventListener('DOMContentLoaded', () => {
    // Add hover effects to cards
    const cards = document.querySelectorAll('.feature-card, .tech-card, .sample-btn');
    cards.forEach(card => {
        card.addEventListener('mouseenter', function() {
            this.style.transform = 'translateY(-2px)';
        });

        card.addEventListener('mouseleave', function() {
            this.style.transform = 'translateY(0)';
        });
    });

    // Add smooth scrolling for internal links
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function (e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));
            if (target) {
                target.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
            }
        });
    });

    // Add loading animation to buttons
    const buttons = document.querySelectorAll('button');
    buttons.forEach(button => {
        button.addEventListener('click', function() {
            if (!this.disabled) {
                this.style.transform = 'scale(0.98)';
                setTimeout(() => {
                    this.style.transform = 'scale(1)';
                }, 100);
            }
        });
    });
});