// MD Simulation Pipeline JavaScript
console.log('Script loading...'); // Debug log

class MDSimulationPipeline {
    constructor() {
        this.currentProtein = null;
        this.preparedProtein = null;
        this.completedProtein = null;
        this.missingResiduesInfo = null;
        this.missingResiduesPdbId = null;
        this.chainSequences = null;
        this.simulationParams = {};
        this.generatedFiles = {};
        this.nglStage = null;
        this.preparedNglStage = null;
        this.completedNglStage = null;
        this.originalNglStage = null;
        this.currentRepresentation = 'cartoon';
        this.preparedRepresentation = 'cartoon';
        this.completedRepresentation = 'cartoon';
        this.originalRepresentation = 'cartoon';
        this.isSpinning = false;
        this.preparedIsSpinning = false;
        this.completedIsSpinning = false;
        this.originalIsSpinning = false;
        this.currentTabIndex = 0;
        this.tabOrder = ['protein-loading', 'fill-missing', 'structure-prep', 'simulation-params', 'simulation-steps', 'file-generation'];
        // Consistent chain color palette - same colors for same chain IDs throughout
        this.chainColorPalette = [
            '#1f77b4', // blue
            '#ff7f0e', // orange
            '#2ca02c', // green
            '#d62728', // red
            '#9467bd', // purple
            '#8c564b', // brown
            '#e377c2', // pink
            '#7f7f7f', // gray
            '#bcbd22', // olive
            '#17becf', // cyan
            '#aec7e8', // light blue
            '#ffbb78', // light orange
            '#98df8a', // light green
            '#ff9896', // light red
            '#c5b0d5', // light purple
            '#c49c94', // light brown
            '#f7b6d3', // light pink
            '#c7c7c7', // light gray
            '#dbdb8d', // light olive
            '#9edae5'  // light cyan
        ];
        this.chainColorMap = {}; // Will store chain ID -> color mapping
        this.init();
        this.initializeTooltips();
    }

    init() {
        this.setupEventListeners();
        this.initializeTabs();
        this.initializeStepToggles();
        this.loadDefaultParams();
        this.updateNavigationState();
    }

    initializeTooltips() {
        // Initialize Bootstrap tooltips using vanilla JavaScript
        // Note: This requires Bootstrap to be loaded
        if (typeof bootstrap !== 'undefined' && bootstrap.Tooltip) {
            const tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-toggle="tooltip"]'));
            tooltipTriggerList.map(function (tooltipTriggerEl) {
                return new bootstrap.Tooltip(tooltipTriggerEl);
            });
        } else {
            console.log('Bootstrap not loaded, tooltips will not work');
        }
    }

    setupEventListeners() {
        // Tab navigation
        document.querySelectorAll('.tab-button').forEach(button => {
            button.addEventListener('click', (e) => this.switchTab(e.target.dataset.tab));
        });

        // File upload
        const fileInput = document.getElementById('pdb-file');
        const fileUploadArea = document.getElementById('file-upload-area');
        const chooseFileBtn = document.getElementById('choose-file-btn');
        
        console.log('File input element:', fileInput);
        console.log('File upload area:', fileUploadArea);
        console.log('Choose file button:', chooseFileBtn);
        
        if (!fileInput) {
            console.error('File input element not found!');
            return;
        }
        
        fileInput.addEventListener('change', (e) => this.handleFileUpload(e));
        
        // Handle click on upload area (but not on the button)
        fileUploadArea.addEventListener('click', (e) => {
            // Only trigger if not clicking on the button
            if (e.target !== chooseFileBtn && !chooseFileBtn.contains(e.target)) {
                console.log('Upload area clicked, triggering file input');
                fileInput.click();
            }
        });
        
        // Handle click on choose file button
        chooseFileBtn.addEventListener('click', (e) => {
            e.stopPropagation(); // Prevent triggering the upload area click
            console.log('Choose file button clicked, triggering file input');
            fileInput.click();
        });
        
        fileUploadArea.addEventListener('dragover', (e) => this.handleDragOver(e));
        fileUploadArea.addEventListener('drop', (e) => this.handleDrop(e));

        // PDB fetch
        document.getElementById('fetch-pdb').addEventListener('click', () => this.fetchPDB());

        // Missing residues analysis
        const detectMissingBtn = document.getElementById('detect-missing-residues');
        if (detectMissingBtn) {
            detectMissingBtn.addEventListener('click', () => this.detectMissingResidues());
        }
        const buildCompleteBtn = document.getElementById('build-complete-structure');
        if (buildCompleteBtn) {
            buildCompleteBtn.addEventListener('click', () => this.buildCompletedStructure());
        }
        const applyTrimBtn = document.getElementById('apply-trim');
        if (applyTrimBtn) {
            applyTrimBtn.addEventListener('click', () => this.applyTrimming());
        }
        const previewCompletedBtn = document.getElementById('preview-completed-structure');
        if (previewCompletedBtn) {
            previewCompletedBtn.addEventListener('click', () => this.previewCompletedStructure());
        }
        const previewSuperimposedBtn = document.getElementById('preview-superimposed-structure');
        if (previewSuperimposedBtn) {
            previewSuperimposedBtn.addEventListener('click', () => this.previewSuperimposedStructure());
        }
        const downloadCompletedBtn = document.getElementById('download-completed-structure');
        if (downloadCompletedBtn) {
            downloadCompletedBtn.addEventListener('click', () => this.downloadCompletedStructure());
        }

        // File generation
        document.getElementById('generate-files').addEventListener('click', () => this.generateAllFiles());
        document.getElementById('preview-files').addEventListener('click', () => this.previewFiles());
        document.getElementById('preview-solvated').addEventListener('click', () => this.previewSolvatedProtein());
        document.getElementById('download-zip').addEventListener('click', () => this.downloadZip());
        

        // Structure preparation
        document.getElementById('prepare-structure').addEventListener('click', () => this.prepareStructure());
        document.getElementById('preview-prepared').addEventListener('click', () => this.previewPreparedStructure());
        document.getElementById('download-prepared').addEventListener('click', () => this.downloadPreparedStructure());
        
        // Ligand download button
        const downloadLigandBtn = document.getElementById('download-ligand');
        if (downloadLigandBtn) {
            downloadLigandBtn.addEventListener('click', (e) => {
                e.preventDefault();
                this.downloadLigandFile();
            });
        }

        // Navigation buttons
        document.getElementById('prev-tab').addEventListener('click', () => this.previousTab());
        document.getElementById('next-tab').addEventListener('click', () => this.nextTab());

        // Parameter changes
        document.querySelectorAll('input, select').forEach(input => {
            input.addEventListener('change', () => this.updateSimulationParams());
        });

        // Render chain and ligand choices when structure tab becomes visible
        document.querySelector('[data-tab="structure-prep"]').addEventListener('click', () => {
            this.renderChainAndLigandSelections();
        });

        // Separate ligands checkbox change
        document.getElementById('separate-ligands').addEventListener('change', (e) => {
            const downloadBtn = document.getElementById('download-ligand');
            
            if (e.target.checked && this.preparedProtein && this.preparedProtein.ligand_present && this.preparedProtein.ligand_content) {
                downloadBtn.disabled = false;
                downloadBtn.classList.remove('btn-outline-secondary');
                downloadBtn.classList.add('btn-outline-primary');
            } else {
                downloadBtn.disabled = true;
                downloadBtn.classList.remove('btn-outline-primary');
                downloadBtn.classList.add('btn-outline-secondary');
            }
        });

        // Preserve ligands checkbox change
        document.getElementById('preserve-ligands').addEventListener('change', (e) => {
            this.toggleLigandForceFieldGroup(e.target.checked);
        });
    }

    initializeTabs() {
        const tabs = document.querySelectorAll('.tab-content');
        tabs.forEach(tab => {
            if (!tab.classList.contains('active')) {
                tab.style.display = 'none';
            }
        });
    }

    initializeStepToggles() {
        document.querySelectorAll('.step-header').forEach(header => {
            header.addEventListener('click', () => {
                const stepItem = header.parentElement;
                const content = stepItem.querySelector('.step-content');
                const isActive = content.classList.contains('active');
                
                // Close all other step contents
                document.querySelectorAll('.step-content').forEach(c => c.classList.remove('active'));
                
                // Toggle current step
                if (!isActive) {
                    content.classList.add('active');
                }
            });
        });
    }

    loadDefaultParams() {
        this.simulationParams = {
            boxType: 'cubic',
            boxSize: 1.0,
            boxMargin: 1.0,
            forceField: 'amber99sb-ildn',
            waterModel: 'tip3p',
            ionConcentration: 150,
            temperature: 300,
            pressure: 1.0,
            couplingType: 'berendsen',
            timestep: 0.002,
            cutoff: 1.0,
            pmeOrder: 4,
            steps: {
                restrainedMin: { enabled: true, steps: 1000, force: 1000 },
                minimization: { enabled: true, steps: 5000, algorithm: 'steep' },
                nvt: { enabled: true, steps: 50000, temperature: 300 },
                npt: { enabled: true, steps: 100000, temperature: 300, pressure: 1.0 },
                production: { enabled: true, steps: 1000000, temperature: 300, pressure: 1.0 }
            }
        };
    }

    switchTab(tabName) {
        // Hide all tab contents
        document.querySelectorAll('.tab-content').forEach(tab => {
            tab.classList.remove('active');
            tab.style.display = 'none';
        });

        // Remove active class from all tab buttons
        document.querySelectorAll('.tab-button').forEach(button => {
            button.classList.remove('active');
        });

        // Show selected tab
        document.getElementById(tabName).classList.add('active');
        document.getElementById(tabName).style.display = 'block';
        
        // Add active class to clicked button
        document.querySelector(`[data-tab="${tabName}"]`).classList.add('active');

        // Update current tab index and navigation state
        this.currentTabIndex = this.tabOrder.indexOf(tabName);
        this.updateNavigationState();
    }

    previousTab() {
        if (this.currentTabIndex > 0) {
            const prevTab = this.tabOrder[this.currentTabIndex - 1];
            this.switchTab(prevTab);
        }
    }

    nextTab() {
        if (this.currentTabIndex < this.tabOrder.length - 1) {
            const nextTab = this.tabOrder[this.currentTabIndex + 1];
            this.switchTab(nextTab);
        }
    }

    updateNavigationState() {
        const prevBtn = document.getElementById('prev-tab');
        const nextBtn = document.getElementById('next-tab');
        const currentStepSpan = document.getElementById('current-step');
        const totalStepsSpan = document.getElementById('total-steps');

        // Update button states
        prevBtn.disabled = this.currentTabIndex === 0;
        nextBtn.disabled = this.currentTabIndex === this.tabOrder.length - 1;

        // Update step indicator
        if (currentStepSpan) {
            currentStepSpan.textContent = this.currentTabIndex + 1;
        }
        if (totalStepsSpan) {
            totalStepsSpan.textContent = this.tabOrder.length;
        }

        // Update next button text based on current tab
        if (this.currentTabIndex === this.tabOrder.length - 1) {
            nextBtn.innerHTML = 'Complete <i class="fas fa-check"></i>';
        } else {
            nextBtn.innerHTML = 'Next <i class="fas fa-chevron-right"></i>';
        }
    }

    handleDragOver(e) {
        e.preventDefault();
        e.currentTarget.style.background = '#e3f2fd';
    }

    handleDrop(e) {
        e.preventDefault();
        e.currentTarget.style.background = '#f8f9fa';
        
        const files = e.dataTransfer.files;
        if (files.length > 0) {
            this.processFile(files[0]);
        }
    }

    handleFileUpload(e) {
        console.log('File upload triggered');
        console.log('Files:', e.target.files);
        const file = e.target.files[0];
        if (file) {
            console.log('File selected:', file.name, file.size, file.type);
            this.processFile(file);
        } else {
            console.log('No file selected');
        }
    }

    processFile(file) {
        console.log('Processing file:', file.name, file.size, file.type);
        
        if (!file.name.toLowerCase().endsWith('.pdb') && !file.name.toLowerCase().endsWith('.ent')) {
            console.log('Invalid file type:', file.name);
            this.showStatus('error', 'Please upload a valid PDB file (.pdb or .ent)');
            return;
        }

        console.log('File validation passed, reading file...');
        const reader = new FileReader();
        reader.onload = (e) => {
            console.log('File read successfully, content length:', e.target.result.length);
            const content = e.target.result;
            this.parsePDBFile(content, file.name);
        };
        reader.onerror = (e) => {
            console.error('Error reading file:', e);
            this.showStatus('error', 'Error reading file');
        };
        reader.readAsText(file);
    }

    async parsePDBFile(content, filename) {
        return this._parsePDBFileInternal(content, filename, true);
    }

    async _parsePDBFileInternal(content, filename, cleanOutput = false) {
        try {
            // Clean output folder when new PDB is loaded (only if requested)
            if (cleanOutput) {
                try {
                    await fetch('/api/clean-output', { method: 'POST' });
                } catch (error) {
                    console.log('Could not clean output folder:', error);
                }
            }
            
            // Save the PDB file to output directory for backend processing
            try {
                await fetch('/api/save-pdb-file', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        pdb_content: content,
                        filename: filename
                    })
                });
            } catch (error) {
                console.log('Could not save PDB file to output directory:', error);
            }

            const lines = content.split('\n');
            let atomCount = 0;
            let chains = new Set();
            let residues = new Set();
            let waterMolecules = 0;
            let ions = 0;
            let ligands = new Set();
            let ligandDetails = [];
            let ligandGroups = new Map(); // Group by RESN-CHAIN
            let hetatoms = 0;
            let structureId = filename.replace(/\.(pdb|ent)$/i, '').toUpperCase();
            
            // Common water molecule names
            const waterNames = new Set(['HOH', 'WAT', 'TIP3', 'TIP4', 'SPC', 'SPCE']);
            
            // Common ion names (expanded to include more ions)
            const ionNames = new Set(['NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'NI', 'CO', 
                'CD', 'HG', 'PB', 'SR', 'BA', 'RB', 'CS', 'LI', 'F', 'BR', 'I', 'SO4', 'PO4', 'CO3', 'NO3', 'NH4']);
            
            // Track ligand entities (separated by TER or different chains)
            let ligandEntities = new Map(); // Map to store ligand entities
            let currentLigandEntity = null;
            let currentChain = null;
            let currentResidue = null;
            
            // Track unique water molecules by residue
            const uniqueWaterResidues = new Set();

            lines.forEach(line => {
                if (line.startsWith('ATOM')) {
                    atomCount++;
                    const chainId = line.substring(21, 22).trim();
                    if (chainId) chains.add(chainId);
                    
                    const resName = line.substring(17, 20).trim();
                    const resNum = line.substring(22, 26).trim();
                    residues.add(`${resName}${resNum}`);
                } else if (line.startsWith('HETATM')) {
                    hetatoms++;
                    const resName = line.substring(17, 20).trim();
                    const resNum = line.substring(22, 26).trim();
                    const chainId = line.substring(21, 22).trim();
                    const entityKey = `${resName}_${resNum}_${chainId}`;
                    
                    if (waterNames.has(resName)) {
                        waterMolecules++;
                        uniqueWaterResidues.add(entityKey);
                    } else if (ionNames.has(resName)) {
                        ions++;
                    } else {
                        // Everything else is treated as ligand (FALLBACK LOGIC)
                        ligands.add(resName);
                        ligandDetails.push({ resn: resName, chain: chainId, resi: resNum });
                        
                        // Group by RESN-CHAIN for UI display
                        const groupKey = `${resName}-${chainId}`;
                        if (!ligandGroups.has(groupKey)) {
                            ligandGroups.set(groupKey, { resn: resName, chain: chainId });
                        }
                        
                        // Track ligand entities (separated by TER or different chains)
                        if (currentChain !== chainId || currentResidue !== resName) {
                            // New ligand entity detected
                            currentLigandEntity = `${resName}_${chainId}`;
                            currentChain = chainId;
                            currentResidue = resName;
                            
                            if (!ligandEntities.has(currentLigandEntity)) {
                                ligandEntities.set(currentLigandEntity, {
                                    name: resName,
                                    chain: chainId,
                                    residueNum: resNum,
                                    atomCount: 0
                                });
                            }
                        }
                        
                        // Increment atom count for this ligand entity
                        if (ligandEntities.has(currentLigandEntity)) {
                            ligandEntities.get(currentLigandEntity).atomCount++;
                        }
                    }
                } else if (line.startsWith('TER')) {
                    // TER record separates entities, reset current ligand tracking
                    currentLigandEntity = null;
                    currentChain = null;
                    currentResidue = null;
                }
            });

            // Count unique water molecules
            const uniqueWaterCount = uniqueWaterResidues.size;
            
            // Count ligand entities (different ligands or same ligand in different chains)
            const ligandEntityCount = ligandEntities.size;
            
            // Get unique ligand names
            const uniqueLigandNames = Array.from(ligands);
            
            // Create ligand info string
            let ligandInfo = 'None';
            if (uniqueLigandNames.length > 0) {
                if (ligandEntityCount > 1) {
                    // Multiple entities - show count and names
                    ligandInfo = `${ligandEntityCount} entities: ${uniqueLigandNames.join(', ')}`;
                } else {
                    // Single entity
                    ligandInfo = uniqueLigandNames.join(', ');
                }
            }

            this.currentProtein = {
                filename: filename,
                structureId: structureId,
                atomCount: atomCount,
                chains: Array.from(chains),
                residueCount: residues.size,
                waterMolecules: uniqueWaterCount,
                ions: ions,
                ligands: uniqueLigandNames,
                ligandDetails: ligandDetails,
                ligandGroups: Array.from(ligandGroups.values()), // RESN-CHAIN groups for UI
                ligandEntities: ligandEntityCount,
                ligandInfo: ligandInfo,
                hetatoms: hetatoms,
                content: content
            };

            this.displayProteinInfo();
            this.showStatus('success', `Successfully loaded ${filename}`);
        } catch (error) {
            this.showStatus('error', 'Error parsing PDB file: ' + error.message);
        }
    }

    displayProteinInfo() {
        if (!this.currentProtein) return;

        document.getElementById('structure-id').textContent = this.currentProtein.structureId;
        document.getElementById('atom-count').textContent = this.currentProtein.atomCount.toLocaleString();
        document.getElementById('chain-info').textContent = this.currentProtein.chains.join(', ');
        document.getElementById('residue-count').textContent = this.currentProtein.residueCount.toLocaleString();
        document.getElementById('water-count').textContent = this.currentProtein.waterMolecules.toLocaleString();
        document.getElementById('ion-count').textContent = this.currentProtein.ions.toLocaleString();
        document.getElementById('ligand-info').textContent = this.currentProtein.ligandInfo;
        document.getElementById('hetatm-count').textContent = this.currentProtein.hetatoms.toLocaleString();

        // Build consistent chain color mapping based on chain IDs from Step 1
        this.buildChainColorMap(this.currentProtein.chains);

        document.getElementById('protein-preview').style.display = 'block';
        
        // Load 3D visualization
        this.load3DVisualization();

        // Also refresh chain/ligand lists when protein info is displayed
        this.renderChainAndLigandSelections();
    }

    buildChainColorMap(chains) {
        // Create a consistent mapping of chain IDs to colors
        // Sort chains to ensure consistent ordering
        const sortedChains = [...chains].sort();
        this.chainColorMap = {};
        sortedChains.forEach((chain, index) => {
            this.chainColorMap[chain] = this.chainColorPalette[index % this.chainColorPalette.length];
        });
        console.log('Chain color map built:', this.chainColorMap);
    }

    async fetchPDB() {
        const pdbId = document.getElementById('pdb-id').value.trim().toUpperCase();
        if (!pdbId) {
            this.showStatus('error', 'Please enter a PDB ID');
            return;
        }

        if (!/^[0-9A-Z]{4}$/.test(pdbId)) {
            this.showStatus('error', 'Please enter a valid 4-character PDB ID');
            return;
        }

        this.showStatus('info', 'Fetching PDB structure...');
        
        try {
            const response = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
            if (!response.ok) {
                throw new Error(`PDB ID ${pdbId} not found`);
            }
            
            const content = await response.text();
            this.parsePDBFile(content, `${pdbId}.pdb`);
            this.showStatus('success', `Successfully fetched PDB structure ${pdbId}`);
        } catch (error) {
            this.showStatus('error', `Error fetching PDB: ${error.message}`);
        }
    }

    showStatus(type, message) {
        const statusDiv = document.getElementById('pdb-status');
        statusDiv.className = `status-message ${type}`;
        statusDiv.textContent = message;
        statusDiv.style.display = 'block';

        // Auto-hide after 5 seconds for success messages
        if (type === 'success') {
            setTimeout(() => {
                statusDiv.style.display = 'none';
            }, 5000);
        }
    }

    updateSimulationParams() {
        // Update basic parameters
        this.simulationParams.boxType = document.getElementById('box-type').value;
        this.simulationParams.boxSize = parseFloat(document.getElementById('box-size').value);
        this.simulationParams.forceField = document.getElementById('force-field').value;
        this.simulationParams.waterModel = document.getElementById('water-model').value;
        this.simulationParams.addIons = document.getElementById('add-ions').value;
        this.simulationParams.temperature = parseInt(document.getElementById('temperature').value);
        this.simulationParams.pressure = parseFloat(document.getElementById('pressure').value);
        this.simulationParams.couplingType = document.getElementById('coupling-type').value;
        this.simulationParams.timestep = parseFloat(document.getElementById('timestep').value);
        this.simulationParams.cutoff = parseFloat(document.getElementById('cutoff').value);
        this.simulationParams.electrostatic = document.getElementById('electrostatic').value;
        this.simulationParams.ligandForceField = document.getElementById('ligand-forcefield').value;

        // Update step parameters
        this.simulationParams.steps.restrainedMin = {
            enabled: document.getElementById('enable-restrained-min').checked,
            steps: parseInt(document.getElementById('restrained-steps').value),
            force: parseInt(document.getElementById('restrained-force').value)
        };

        this.simulationParams.steps.minimization = {
            enabled: document.getElementById('enable-minimization').checked,
            steps: parseInt(document.getElementById('min-steps').value),
            algorithm: document.getElementById('min-algorithm').value
        };

        this.simulationParams.steps.nvt = {
            enabled: document.getElementById('enable-nvt').checked,
            steps: parseInt(document.getElementById('nvt-steps').value),
            temperature: parseInt(document.getElementById('nvt-temp').value)
        };

        this.simulationParams.steps.npt = {
            enabled: document.getElementById('enable-npt').checked,
            steps: parseInt(document.getElementById('npt-steps').value),
            temperature: parseInt(document.getElementById('npt-temp').value),
            pressure: parseFloat(document.getElementById('npt-pressure').value)
        };

        this.simulationParams.steps.production = {
            enabled: document.getElementById('enable-production').checked,
            steps: parseInt(document.getElementById('prod-steps').value),
            temperature: parseInt(document.getElementById('prod-temp').value),
            pressure: parseFloat(document.getElementById('prod-pressure').value)
        };
    }

    toggleLigandForceFieldGroup(show) {
        const section = document.getElementById('ligand-forcefield-section');
        if (show) {
            section.style.display = 'block';
            section.classList.remove('disabled');
        } else {
            section.style.display = 'none';
            section.classList.add('disabled');
        }
    }

    async calculateNetCharge(event) {
        console.log('calculateNetCharge called'); // Debug log
        if (!this.preparedProtein) {
            alert('Please prepare structure first before calculating net charge.');
            return;
        }

        // Show loading state
        const button = event ? event.target : document.querySelector('button[onclick*="calculateNetCharge"]');
        const originalText = button.innerHTML;
        button.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Calculating...';
        button.disabled = true;

        try {
            // Get the selected force field
            const selectedForceField = document.getElementById('force-field').value;
            
            const response = await fetch('/api/calculate-net-charge', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    force_field: selectedForceField
                })
            });

            const result = await response.json();

            if (result.success) {
                // Update the Add Ions dropdown based on suggestion
                const addIonsSelect = document.getElementById('add-ions');
                if (result.ion_type === 'Cl-') {
                    addIonsSelect.value = 'Cl-';
                } else if (result.ion_type === 'Na+') {
                    addIonsSelect.value = 'Na+';
                } else {
                    addIonsSelect.value = 'None';
                }

                // Show detailed results
                alert(`✅ System Charge Analysis Complete!\n\n` +
                      `Net Charge: ${result.net_charge}\n` +
                      `Recommendation: ${result.suggestion}\n` +
                      `Ligand Present: ${result.ligand_present ? 'Yes' : 'No'}`);
            } else {
                alert(`❌ Error: ${result.error}`);
            }
        } catch (error) {
            console.error('Error calculating net charge:', error);
            alert(`❌ Error: Failed to calculate net charge. ${error.message}`);
        } finally {
            // Restore button state
            button.innerHTML = originalText;
            button.disabled = false;
        }
    }

    async generateLigandFF(event) {
        console.log('generateLigandFF called'); // Debug log
        if (!this.preparedProtein || !this.preparedProtein.ligand_present) {
            alert('No ligand found. Please ensure ligands are preserved during structure preparation.');
            return;
        }

        const selectedFF = document.getElementById('ligand-forcefield').value;
        
        // Show loading state
        const button = event ? event.target : document.querySelector('button[onclick*="generateLigandFF"]');
        const originalText = button.innerHTML;
        button.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Generating...';
        button.disabled = true;

        try {
            const response = await fetch('/api/generate-ligand-ff', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    force_field: selectedFF
                })
            });

            const result = await response.json();

            if (result.success) {
                alert(`✅ ${result.message}\n\nNet charge: ${result.net_charge}\n\nGenerated files:\n- ${result.files.mol2}\n- ${result.files.frcmod}`);
            } else {
                alert(`❌ Error: ${result.error}`);
            }
        } catch (error) {
            console.error('Error generating ligand force field:', error);
            alert(`❌ Error: Failed to generate force field parameters. ${error.message}`);
        } finally {
            // Restore button state
            button.innerHTML = originalText;
            button.disabled = false;
        }
    }

    countAtomsInPDB(pdbContent) {
        const lines = pdbContent.split('\n');
        return lines.filter(line => line.startsWith('ATOM') || line.startsWith('HETATM')).length;
    }

    async generateAllFiles() {
        if (!this.preparedProtein) {
            alert('Please prepare structure first');
            return;
        }

        // Show loading state
        const button = document.getElementById('generate-files');
        const originalText = button.innerHTML;
        button.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Generating...';
        button.disabled = true;

        try {
               // Collect all simulation parameters
               const params = {
                   cutoff_distance: parseFloat(document.getElementById('cutoff').value),
                   temperature: parseFloat(document.getElementById('temperature').value),
                   pressure: parseFloat(document.getElementById('pressure').value),
                   restrained_steps: parseInt(document.getElementById('restrained-steps').value),
                   restrained_force: parseFloat(document.getElementById('restrained-force').value),
                   min_steps: parseInt(document.getElementById('min-steps').value),
                   npt_heating_steps: parseInt(document.getElementById('nvt-steps').value),
                   npt_equilibration_steps: parseInt(document.getElementById('npt-steps').value),
                   production_steps: parseInt(document.getElementById('prod-steps').value),
                   timestep: parseFloat(document.getElementById('timestep').value),
                   // Force field parameters
                   force_field: document.getElementById('force-field').value,
                   water_model: document.getElementById('water-model').value,
                   add_ions: document.getElementById('add-ions').value,
                   distance: parseFloat(document.getElementById('box-size').value)
               };

            const response = await fetch('/api/generate-all-files', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(params)
            });

            const result = await response.json();

            if (result.success) {
                let message = `✅ ${result.message}\n\nGenerated files:\n`;
                result.files_generated.forEach(file => {
                    message += `- ${file}\n`;
                });
                
                if (result.warnings && result.warnings.length > 0) {
                    message += `\n⚠️ Warnings:\n`;
                    result.warnings.forEach(warning => {
                        message += `- ${warning}\n`;
                    });
                }
                
                alert(message);

                // Reveal the download section
                const downloadSection = document.getElementById('download-section');
                if (downloadSection) {
                    downloadSection.style.display = 'block';
                }
            } else {
                alert(`❌ Error: ${result.error}`);
            }
        } catch (error) {
            console.error('Error generating files:', error);
            alert(`❌ Error: Failed to generate simulation files. ${error.message}`);
        } finally {
            // Restore button state
            button.innerHTML = originalText;
            button.disabled = false;
        }
    }

    createSimulationFiles() {
        const files = {};
        const proteinName = this.currentProtein.structureId.toLowerCase();
        
        // Generate GROMACS input files
        files[`${proteinName}.mdp`] = this.generateMDPFile();
        files[`${proteinName}_restrained.mdp`] = this.generateRestrainedMDPFile();
        files[`${proteinName}_min.mdp`] = this.generateMinimizationMDPFile();
        files[`${proteinName}_nvt.mdp`] = this.generateNVTMDPFile();
        files[`${proteinName}_npt.mdp`] = this.generateNPTMDPFile();
        files[`${proteinName}_prod.mdp`] = this.generateProductionMDPFile();
        
        // Generate PBS script
        files[`${proteinName}_simulation.pbs`] = this.generatePBSScript();
        
        // Generate setup script
        files[`setup_${proteinName}.sh`] = this.generateSetupScript();
        
        // Generate analysis script
        files[`analyze_${proteinName}.sh`] = this.generateAnalysisScript();

        return files;
    }

    generateMDPFile() {
        const params = this.simulationParams;
        return `; MD Simulation Parameters
; Generated by MD Simulation Pipeline

; Run parameters
integrator = md
dt = ${params.timestep}
nsteps = ${params.steps.production.steps}

; Output control
nstxout = 5000
nstvout = 5000
nstenergy = 1000
nstlog = 1000

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter = 1
lincs_order = 4

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 40
rlist = ${params.cutoff}

; Electrostatics
coulombtype = PME
rcoulomb = ${params.cutoff}
pme_order = ${params.pmeOrder}
fourierspacing = 0.16

; Van der Waals
vdwtype = Cut-off
rvdw = ${params.cutoff}

; Temperature coupling
tcoupl = ${params.couplingType}
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = ${params.temperature} ${params.temperature}

; Pressure coupling
pcoupl = ${params.couplingType}
pcoupltype = isotropic
tau_p = 2.0
ref_p = ${params.pressure}
compressibility = 4.5e-5

; Dispersion correction
DispCorr = EnerPres

; Velocity generation
gen_vel = yes
gen_temp = ${params.temperature}
gen_seed = -1
`;
    }

    generateRestrainedMDPFile() {
        const params = this.simulationParams;
        return `; Restrained Minimization Parameters
integrator = steep
nsteps = ${params.steps.restrainedMin.steps}
emstep = 0.01
emtol = 1000

; Position restraints
define = -DPOSRES
refcoord_scaling = com

; Output control
nstxout = 100
nstenergy = 100
nstlog = 100

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = ${params.cutoff}

; Electrostatics
coulombtype = PME
rcoulomb = ${params.cutoff}
pme_order = ${params.pme_order}

; Van der Waals
vdwtype = Cut-off
rvdw = ${params.cutoff}
`;
    }

    generateMinimizationMDPFile() {
        const params = this.simulationParams;
        return `; Minimization Parameters
integrator = ${params.steps.minimization.algorithm}
nsteps = ${params.steps.minimization.steps}
emstep = 0.01
emtol = 1000

; Output control
nstxout = 100
nstenergy = 100
nstlog = 100

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 10
rlist = ${params.cutoff}

; Electrostatics
coulombtype = PME
rcoulomb = ${params.cutoff}
pme_order = ${params.pme_order}

; Van der Waals
vdwtype = Cut-off
rvdw = ${params.cutoff}
`;
    }

    generateNVTMDPFile() {
        const params = this.simulationParams;
        return `; NVT Equilibration Parameters
integrator = md
dt = ${params.timestep}
nsteps = ${params.steps.nvt.steps}

; Output control
nstxout = 5000
nstvout = 5000
nstenergy = 1000
nstlog = 1000

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter = 1
lincs_order = 4

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 40
rlist = ${params.cutoff}

; Electrostatics
coulombtype = PME
rcoulomb = ${params.cutoff}
pme_order = ${params.pme_order}

; Van der Waals
vdwtype = Cut-off
rvdw = ${params.cutoff}

; Temperature coupling
tcoupl = ${params.couplingType}
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = ${params.steps.nvt.temperature} ${params.steps.nvt.temperature}

; Pressure coupling (disabled for NVT)
pcoupl = no

; Velocity generation
gen_vel = yes
gen_temp = ${params.steps.nvt.temperature}
gen_seed = -1
`;
    }

    generateNPTMDPFile() {
        const params = this.simulationParams;
        return `; NPT Equilibration Parameters
integrator = md
dt = ${params.timestep}
nsteps = ${params.steps.npt.steps}

; Output control
nstxout = 5000
nstvout = 5000
nstenergy = 1000
nstlog = 1000

; Bond parameters
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter = 1
lincs_order = 4

; Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
nstlist = 40
rlist = ${params.cutoff}

; Electrostatics
coulombtype = PME
rcoulomb = ${params.cutoff}
pme_order = ${params.pme_order}

; Van der Waals
vdwtype = Cut-off
rvdw = ${params.cutoff}

; Temperature coupling
tcoupl = ${params.couplingType}
tc-grps = Protein Non-Protein
tau_t = 0.1 0.1
ref_t = ${params.steps.npt.temperature} ${params.steps.npt.temperature}

; Pressure coupling
pcoupl = ${params.couplingType}
pcoupltype = isotropic
tau_p = 2.0
ref_p = ${params.steps.npt.pressure}
compressibility = 4.5e-5

; Velocity generation
gen_vel = no
`;
    }

    generateProductionMDPFile() {
        return this.generateMDPFile(); // Same as main MDP file
    }

    generatePBSScript() {
        const proteinName = this.currentProtein.structureId.toLowerCase();
        const totalSteps = this.simulationParams.steps.production.steps;
        const timeInNs = (totalSteps * this.simulationParams.timestep) / 1000;
        
        return `#!/bin/bash
#PBS -N ${proteinName}_md
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -q normal
#PBS -j oe

# Change to the directory where the job was submitted
cd $PBS_O_WORKDIR

# Load required modules
module load gromacs/2023.2
module load intel/2021.4.0

# Set up environment
export OMP_NUM_THREADS=16
export GMX_MAXBACKUP=-1

# Simulation parameters
PROTEIN=${proteinName}
STEPS=${totalSteps}
TIME_NS=${timeInNs.toFixed(2)}

echo "Starting MD simulation for $PROTEIN"
echo "Total simulation time: $TIME_NS ns"
echo "Job started at: $(date)"

# Run the simulation
./run_simulation.sh $PROTEIN

echo "Simulation completed at: $(date)"
echo "Results saved in output directory"
`;
    }

    generateSetupScript() {
        const proteinName = this.currentProtein.structureId.toLowerCase();
        return `#!/bin/bash
# Setup script for ${proteinName} MD simulation
# Generated by MD Simulation Pipeline

set -e

PROTEIN=${proteinName}
FORCE_FIELD=${this.simulationParams.forceField}
WATER_MODEL=${this.simulationParams.waterModel}

echo "Setting up MD simulation for $PROTEIN"

# Create output directory
mkdir -p output

# 1. Prepare protein structure
echo "Preparing protein structure..."
gmx pdb2gmx -f ${PROTEIN}.pdb -o ${PROTEIN}_processed.gro -p ${PROTEIN}.top -ff ${FORCE_FIELD} -water ${WATER_MODEL}

# 2. Define simulation box
echo "Defining simulation box..."
gmx editconf -f ${PROTEIN}_processed.gro -o ${PROTEIN}_box.gro -c -d ${this.simulationParams.boxMargin} -bt ${this.simulationParams.boxType}

# 3. Add solvent
echo "Adding solvent..."
gmx solvate -cp ${PROTEIN}_box.gro -cs spc216.gro -o ${PROTEIN}_solv.gro -p ${PROTEIN}.top

# 4. Add ions
echo "Adding ions..."
gmx grompp -f ${PROTEIN}_restrained.mdp -c ${PROTEIN}_solv.gro -p ${PROTEIN}.top -o ${PROTEIN}_ions.tpr
echo "SOL" | gmx genion -s ${PROTEIN}_ions.tpr -o ${PROTEIN}_final.gro -p ${PROTEIN}.top -pname NA -nname CL -neutral

echo "Setup completed successfully!"
echo "Ready to run simulation with: ./run_simulation.sh $PROTEIN"
`;
    }

    generateAnalysisScript() {
        const proteinName = this.currentProtein.structureId.toLowerCase();
        return `#!/bin/bash
# Analysis script for ${proteinName} MD simulation
# Generated by MD Simulation Pipeline

PROTEIN=${proteinName}

echo "Analyzing MD simulation results for $PROTEIN"

# Create analysis directory
mkdir -p analysis

# 1. RMSD analysis
echo "Calculating RMSD..."
echo "Protein" | gmx rms -s ${PROTEIN}_final.tpr -f ${PROTEIN}_prod.xtc -o analysis/${PROTEIN}_rmsd.xvg -tu ns

# 2. RMSF analysis
echo "Calculating RMSF..."
echo "Protein" | gmx rmsf -s ${PROTEIN}_final.tpr -f ${PROTEIN}_prod.xtc -o analysis/${PROTEIN}_rmsf.xvg -res

# 3. Radius of gyration
echo "Calculating radius of gyration..."
echo "Protein" | gmx gyrate -s ${PROTEIN}_final.tpr -f ${PROTEIN}_prod.xtc -o analysis/${PROTEIN}_gyrate.xvg

# 4. Hydrogen bonds
echo "Analyzing hydrogen bonds..."
echo "Protein" | gmx hbond -s ${PROTEIN}_final.tpr -f ${PROTEIN}_prod.xtc -num analysis/${PROTEIN}_hbonds.xvg

# 5. Energy analysis
echo "Analyzing energies..."
gmx energy -f ${PROTEIN}_prod.edr -o analysis/${PROTEIN}_energy.xvg

# 6. Generate plots
echo "Generating analysis plots..."
python3 plot_analysis.py ${PROTEIN}

echo "Analysis completed! Results saved in analysis/ directory"
`;
    }

    displayGeneratedFiles() {
        const filesList = document.getElementById('files-list');
        filesList.innerHTML = '';

        Object.entries(this.generatedFiles).forEach(([filename, content]) => {
            const fileItem = document.createElement('div');
            fileItem.className = 'file-item';
            
            const fileType = this.getFileType(filename);
            const fileSize = this.formatFileSize(content.length);
            
            fileItem.innerHTML = `
                <h4><i class="fas ${this.getFileIcon(filename)}"></i> ${filename}</h4>
                <p><strong>Type:</strong> ${fileType}</p>
                <p><strong>Size:</strong> ${fileSize}</p>
                <button class="btn btn-secondary btn-sm" onclick="mdPipeline.previewFile('${filename}')">
                    <i class="fas fa-eye"></i> Preview
                </button>
                <button class="btn btn-primary btn-sm" onclick="mdPipeline.downloadFile('${filename}')">
                    <i class="fas fa-download"></i> Download
                </button>
            `;
            
            filesList.appendChild(fileItem);
        });
    }

    getFileType(filename) {
        const extension = filename.split('.').pop().toLowerCase();
        const types = {
            'mdp': 'GROMACS MDP',
            'pbs': 'PBS Script',
            'sh': 'Shell Script',
            'gro': 'GROMACS Structure',
            'top': 'GROMACS Topology',
            'xvg': 'GROMACS Data'
        };
        return types[extension] || 'Text File';
    }

    getFileIcon(filename) {
        const extension = filename.split('.').pop().toLowerCase();
        const icons = {
            'mdp': 'fa-cogs',
            'pbs': 'fa-tasks',
            'sh': 'fa-terminal',
            'gro': 'fa-cube',
            'top': 'fa-sitemap',
            'xvg': 'fa-chart-line'
        };
        return icons[extension] || 'fa-file';
    }

    formatFileSize(bytes) {
        if (bytes < 1024) return bytes + ' B';
        if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
        return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
    }

    previewFile(filename) {
        const content = this.generatedFiles[filename];
        const previewWindow = window.open('', '_blank', 'width=800,height=600');
        previewWindow.document.write(`
            <html>
                <head>
                    <title>Preview: ${filename}</title>
                    <style>
                        body { font-family: monospace; margin: 20px; background: #f5f5f5; }
                        pre { background: white; padding: 20px; border-radius: 5px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }
                        h1 { color: #333; }
                    </style>
                </head>
                <body>
                    <h1>${filename}</h1>
                    <pre>${content}</pre>
                </body>
            </html>
        `);
    }

    downloadFile(filename) {
        const content = this.generatedFiles[filename];
        const blob = new Blob([content], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    async previewFiles() {
        try {
            const resp = await fetch('/api/get-generated-files');
            const data = await resp.json();
            if (!data.success) {
                alert('❌ Error: ' + (data.error || 'Unable to load files'));
                return;
            }
            const filesList = document.getElementById('files-list');
            if (!filesList) return;
            filesList.innerHTML = '';
            
            // Store file contents for modal display
            this.fileContents = data.files;
            
            Object.entries(data.files).forEach(([name, content]) => {
                const fileItem = document.createElement('div');
                fileItem.className = 'file-item';
                fileItem.style.cssText = 'padding: 10px; margin: 5px 0; border: 1px solid #ddd; border-radius: 5px; cursor: pointer; background: #f9f9f9;';
                fileItem.innerHTML = `<strong>${name}</strong>`;
                fileItem.onclick = () => this.showFileContent(name, content);
                filesList.appendChild(fileItem);
            });
            
            // Reveal preview and download areas
            const preview = document.getElementById('files-preview');
            if (preview) preview.style.display = 'block';
            const dl = document.getElementById('download-section');
            if (dl) dl.style.display = 'block';
            this.switchTab('file-generation');
        } catch (e) {
            console.error('Preview error:', e);
            alert('❌ Failed to preview files: ' + e.message);
        }
    }

    showFileContent(filename, content) {
        // Create modal if it doesn't exist
        let modal = document.getElementById('file-content-modal');
        if (!modal) {
            modal = document.createElement('div');
            modal.id = 'file-content-modal';
            modal.style.cssText = `
                position: fixed; top: 0; left: 0; width: 100%; height: 100%; 
                background: rgba(0,0,0,0.5); z-index: 1000; display: none;
            `;
            modal.innerHTML = `
                <div style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);
                           background: white; border-radius: 10px; padding: 20px; max-width: 80%; max-height: 80%;
                           overflow: auto; box-shadow: 0 4px 20px rgba(0,0,0,0.3);">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;">
                        <h3 id="modal-filename" style="margin: 0; color: #333;"></h3>
                        <button id="close-modal" style="background: #dc3545; color: white; border: none; 
                                border-radius: 5px; padding: 8px 15px; cursor: pointer;">Close</button>
                    </div>
                    <pre id="modal-content" style="background: #f8f9fa; padding: 15px; border-radius: 5px; 
                         overflow: auto; max-height: 60vh; white-space: pre-wrap; font-family: monospace;"></pre>
                </div>
            `;
            document.body.appendChild(modal);
            
            // Close modal handlers
            document.getElementById('close-modal').onclick = () => modal.style.display = 'none';
            modal.onclick = (e) => {
                if (e.target === modal) modal.style.display = 'none';
            };
        }
        
        // Populate and show modal
        document.getElementById('modal-filename').textContent = filename;
        document.getElementById('modal-content').textContent = content;
        modal.style.display = 'block';
    }

    async downloadZip() {
        try {
            const resp = await fetch('/api/download-output-zip');
            if (!resp.ok) {
                const text = await resp.text();
                throw new Error(text || 'Failed to create ZIP');
            }
            const blob = await resp.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'output.zip';
            document.body.appendChild(a);
            a.click();
            a.remove();
            window.URL.revokeObjectURL(url);
        } catch (e) {
            console.error('Download error:', e);
            alert('❌ Failed to download ZIP: ' + e.message);
        }
    }

    async previewSolvatedProtein() {
        try {
            // Show loading state
            const button = document.getElementById('preview-solvated');
            const originalText = button.innerHTML;
            button.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Loading...';
            button.disabled = true;

            // Fetch a single viewer PDB that marks ligands as HETATM within protein_solvated frame
            const response = await fetch('/api/get-viewer-pdb');
            if (!response.ok) {
                throw new Error('Viewer PDB not available. Please generate files first.');
            }

            const data = await response.json();
            if (!data.success) {
                throw new Error(data.error || 'Failed to load viewer PDB');
            }

            // Open the dedicated viewer page (bypasses CSP issues)
            window.open('/viewer/viewer_protein_with_ligand.pdb', '_blank');

        } catch (error) {
            console.error('Error previewing solvated protein:', error);
            alert('❌ Error: ' + error.message);
        } finally {
            // Restore button state
            const button = document.getElementById('preview-solvated');
            button.innerHTML = '<i class="fas fa-tint"></i> Preview Solvated Protein';
            button.disabled = false;
        }
    }

    

    displaySimulationSummary() {
        const summaryContent = document.getElementById('summary-content');
        const params = this.simulationParams;
        const protein = this.currentProtein;
        
        const totalTime = (params.steps.production.steps * params.timestep) / 1000; // Convert to ns
        
        summaryContent.innerHTML = `
            <div class="summary-item">
                <h4>Protein Information</h4>
                <p><strong>Structure ID:</strong> ${protein.structureId}</p>
                <p><strong>Atoms:</strong> ${protein.atomCount.toLocaleString()}</p>
                <p><strong>Chains:</strong> ${protein.chains.join(', ')}</p>
                <p><strong>Residues:</strong> ${protein.residueCount.toLocaleString()}</p>
            </div>
            <div class="summary-item">
                <h4>System Components</h4>
                <p><strong>Water molecules:</strong> ${protein.waterMolecules.toLocaleString()}</p>
                <p><strong>Ions:</strong> ${protein.ions.toLocaleString()}</p>
                <p><strong>Ligands:</strong> ${protein.ligands.length > 0 ? protein.ligands.join(', ') : 'None'}</p>
                <p><strong>HETATM entries:</strong> ${protein.hetatoms.toLocaleString()}</p>
            </div>
            <div class="summary-item">
                <h4>Simulation Box</h4>
                <p><strong>Type:</strong> ${params.boxType}</p>
                <p><strong>Size:</strong> ${params.boxSize} nm</p>
                <p><strong>Margin:</strong> ${params.boxMargin} nm</p>
            </div>
            <div class="summary-item">
                <h4>Force Field & Water</h4>
                <p><strong>Force Field:</strong> ${params.forceField}</p>
                <p><strong>Water Model:</strong> ${params.waterModel}</p>
                <p><strong>Ion Conc.:</strong> ${params.ionConcentration} mM</p>
            </div>
            <div class="summary-item">
                <h4>Simulation Parameters</h4>
                <p><strong>Temperature:</strong> ${params.temperature} K</p>
                <p><strong>Pressure:</strong> ${params.pressure} bar</p>
                <p><strong>Time Step:</strong> ${params.timestep} ps</p>
            </div>
            <div class="summary-item">
                <h4>Simulation Time</h4>
                <p><strong>Total Time:</strong> ${totalTime.toFixed(2)} ns</p>
                <p><strong>Steps:</strong> ${params.steps.production.steps.toLocaleString()}</p>
                <p><strong>Output Freq:</strong> Every 5 ps</p>
            </div>
            <div class="summary-item">
                <h4>Generated Files</h4>
                <p><strong>MDP Files:</strong> 6</p>
                <p><strong>Scripts:</strong> 3</p>
                <p><strong>Total Size:</strong> ${this.formatFileSize(Object.values(this.generatedFiles).join('').length)}</p>
            </div>
        `;
    }

    // 3D Visualization Methods
    async load3DVisualization() {
        if (!this.currentProtein) return;

        try {
            // Initialize NGL stage if not already done
            if (!this.nglStage) {
                this.nglStage = new NGL.Stage("ngl-viewer", {
                    backgroundColor: "white",
                    quality: "medium"
                });
            }

            // Clear existing components
            this.nglStage.removeAllComponents();

            // Create a blob from PDB content
            const blob = new Blob([this.currentProtein.content], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);

            // Load the structure
            const component = await this.nglStage.loadFile(url, {
                ext: "pdb",
                defaultRepresentation: false
            });

            // Add cartoon representation for each chain with consistent colors
            // This ensures each chain gets the same color as in Step 1
            // Use chains from parsed protein data (more reliable than structure API)
            const structureChains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                // Add representation for each chain that exists in the structure
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        component.addRepresentation("cartoon", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.9
                        });
                    }
                });
            } else {
                // Fallback: use chainid if color map not available
                component.addRepresentation("cartoon", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.9
                });
            }
            
            // Apply consistent chain colors after representation is added (backup)
            setTimeout(() => {
                this.applyConsistentChainColors(component);
            }, 500);

            // Add ball and stick for water molecules
            if (this.currentProtein.waterMolecules > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "water",
                    color: "cyan",
                    colorScheme: "uniform",
                    radius: 0.1
                });
            }

            // Add ball and stick for ions
            if (this.currentProtein.ions > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "ion",
                    color: "element",
                    radius: 0.2
                });
            }

            // Add ball and stick for ligands
            if (this.currentProtein.ligands.length > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "hetero",
                    color: "element",
                    radius: 0.15
                });
            }

            // Auto-fit the view
            this.nglStage.autoView();

            // Show controls
            document.getElementById('viewer-controls').style.display = 'flex';

            // Clean up the blob URL
            URL.revokeObjectURL(url);

        } catch (error) {
            console.error('Error loading 3D visualization:', error);
            this.showStatus('error', 'Error loading 3D visualization: ' + error.message);
        }
    }

    resetView() {
        if (this.nglStage) {
            this.nglStage.autoView();
        }
    }

    toggleRepresentation() {
        if (!this.nglStage) return;

        const components = this.nglStage.compList;
        if (components.length === 0) return;

        const component = components[0];
        component.removeAllRepresentations();

        if (this.currentRepresentation === 'cartoon') {
            // Switch to ball and stick for everything
            component.addRepresentation("ball+stick", {
                color: "element",
                radius: 0.15
            });
            this.currentRepresentation = 'ball+stick';
            document.getElementById('style-text').textContent = 'Ball & Stick';
        } else if (this.currentRepresentation === 'ball+stick') {
            // Switch to surface (protein only) + ball&stick for others
            // Use consistent chain colors
            // Use chains from parsed protein data (more reliable than structure API)
            const structureChains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        component.addRepresentation("surface", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.7
                        });
                    }
                });
            } else {
                component.addRepresentation("surface", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.7
                });
            }

            // Add ball and stick for water molecules
            if (this.currentProtein.waterMolecules > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "water",
                    color: "cyan",
                    colorScheme: "uniform",
                    radius: 0.1
                });
            }

            // Add ball and stick for ions
            if (this.currentProtein.ions > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "ion",
                    color: "element",
                    radius: 0.2
                });
            }

            // Add ball and stick for ligands
            if (this.currentProtein.ligands.length > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "hetero",
                    color: "element",
                    radius: 0.15
                });
            }

            this.currentRepresentation = 'surface';
            document.getElementById('style-text').textContent = 'Surface';
        } else {
            // Switch back to mixed representation (protein ribbon + others ball&stick)
            component.addRepresentation("cartoon", {
                sele: "protein",
                colorScheme: "chainname",
                opacity: 0.8
            });

            // Add ball and stick for water molecules
            if (this.currentProtein.waterMolecules > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "water",
                    color: "cyan",
                    colorScheme: "uniform",
                    radius: 0.1
                });
            }

            // Add ball and stick for ions
            if (this.currentProtein.ions > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "ion",
                    color: "element",
                    radius: 0.2
                });
            }

            // Add ball and stick for ligands
            if (this.currentProtein.ligands.length > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "hetero",
                    color: "element",
                    radius: 0.15
                });
            }

            this.currentRepresentation = 'cartoon';
            document.getElementById('style-text').textContent = 'Mixed View';
        }
    }

    toggleSpin() {
        if (!this.nglStage) return;

        this.isSpinning = !this.isSpinning;
        this.nglStage.setSpin(this.isSpinning);
    }

    // Structure Preparation Methods
    async prepareStructure() {
        if (!this.currentProtein) {
            alert('Please load a protein structure first');
            return;
        }

        // Get preparation options
        const options = {
            remove_water: document.getElementById('remove-water').checked,
            remove_ions: document.getElementById('remove-ions').checked,
            remove_hydrogens: document.getElementById('remove-hydrogens').checked,
            add_nme: document.getElementById('add-nme').checked,
            add_ace: document.getElementById('add-ace').checked,
            preserve_ligands: document.getElementById('preserve-ligands').checked,
            separate_ligands: document.getElementById('separate-ligands').checked,
            selected_chains: this.getSelectedChains(),
            selected_ligands: this.getSelectedLigands()
        };

        // Show status
        document.getElementById('prep-status').style.display = 'block';
        document.getElementById('prep-status-content').innerHTML = `
            <p><i class="fas fa-spinner fa-spin"></i> Preparing structure...</p>
        `;

        try {
            // Call Python backend
            const response = await fetch('/api/prepare-structure', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    pdb_content: this.currentProtein.content,
                    options: options
                })
            });

            const result = await response.json();

            if (result.success) {
                // Store prepared structure
                this.preparedProtein = {
                    content: result.prepared_structure,
                    original_atoms: result.original_atoms,
                    prepared_atoms: result.prepared_atoms,
                    removed_components: result.removed_components,
                    added_capping: result.added_capping,
                    preserved_ligands: result.preserved_ligands,
                    ligand_present: result.ligand_present,
                    separate_ligands: result.separate_ligands,
                    ligand_content: result.ligand_content || ''
                };

                // Format removed components
                const removedText = result.removed_components ? 
                    Object.entries(result.removed_components)
                        .filter(([key, value]) => value > 0)
                        .map(([key, value]) => `${key}: ${value}`)
                        .join(', ') || 'None' : 'None';
                
                // Format added capping
                const addedText = result.added_capping ? 
                    Object.entries(result.added_capping)
                        .filter(([key, value]) => value > 0)
                        .map(([key, value]) => `${key}: ${value}`)
                        .join(', ') || 'None' : 'None';

                // Update status
                document.getElementById('prep-status-content').innerHTML = `
                    <p><i class="fas fa-check-circle"></i> Structure preparation completed!</p>
                    <p><strong>Original atoms:</strong> ${result.original_atoms.toLocaleString()}</p>
                    <p><strong>Prepared atoms:</strong> ${result.prepared_atoms.toLocaleString()}</p>
                    <p><strong>Removed:</strong> ${removedText}</p>
                    <p><strong>Added:</strong> ${addedText}</p>
                    <p><strong>Ligands:</strong> ${result.preserved_ligands}</p>
                    <p>Ready for AMBER force field generation!</p>
                `;

                // Enable preview and download buttons
                document.getElementById('preview-prepared').disabled = false;
                document.getElementById('download-prepared').disabled = false;
                
                // Enable ligand download button if ligands are present and separate ligands is checked
                const separateLigandsChecked = document.getElementById('separate-ligands').checked;
                const downloadLigandBtn = document.getElementById('download-ligand');
                if (result.ligand_present && separateLigandsChecked && result.ligand_content) {
                    downloadLigandBtn.disabled = false;
                    downloadLigandBtn.classList.remove('btn-outline-secondary');
                    downloadLigandBtn.classList.add('btn-outline-primary');
                } else {
                    downloadLigandBtn.disabled = true;
                    downloadLigandBtn.classList.remove('btn-outline-primary');
                    downloadLigandBtn.classList.add('btn-outline-secondary');
                }

                // Show ligand force field group if preserve ligands is checked
                const preserveLigandsChecked = document.getElementById('preserve-ligands').checked;
                if (preserveLigandsChecked && result.ligand_present) {
                    this.toggleLigandForceFieldGroup(true);
                }
            } else {
                throw new Error(result.error || 'Structure preparation failed');
            }
        } catch (error) {
            console.error('Error preparing structure:', error);
            document.getElementById('prep-status-content').innerHTML = `
                <p><i class="fas fa-exclamation-triangle"></i> Error preparing structure</p>
                <p>${error.message}</p>
            `;
        }
    }

    renderChainAndLigandSelections() {
        if (!this.currentProtein) return;
        // Render chains
        const chainContainer = document.getElementById('chain-selection');
        if (chainContainer) {
            chainContainer.innerHTML = '';
            this.currentProtein.chains.forEach(chainId => {
                const id = `chain-${chainId}`;
                const wrapper = document.createElement('div');
                wrapper.className = 'checkbox-inline';
                wrapper.innerHTML = `
                    <label class="checkbox-container">
                        <input type="checkbox" id="${id}" data-chain="${chainId}">
                        <span class="checkmark"></span>
                        Chain ${chainId}
                    </label>`;
                chainContainer.appendChild(wrapper);
            });
        }

        // Render ligands (RESN-CHAIN groups)
        const ligandContainer = document.getElementById('ligand-selection');
        if (ligandContainer) {
            ligandContainer.innerHTML = '';
            if (Array.isArray(this.currentProtein.ligandGroups) && this.currentProtein.ligandGroups.length > 0) {
                this.currentProtein.ligandGroups.forEach(l => {
                    const key = `${l.resn}-${l.chain}`;
                    const id = `lig-${key}`;
                    const wrapper = document.createElement('div');
                    wrapper.className = 'checkbox-inline';
                    wrapper.innerHTML = `
                        <label class="checkbox-container">
                            <input type="checkbox" id="${id}" data-resn="${l.resn}" data-chain="${l.chain}">
                            <span class="checkmark"></span>
                            ${key}
                        </label>`;
                    ligandContainer.appendChild(wrapper);
                });
            } else {
                // Fallback: show unique ligand names if detailed positions not parsed
                if (Array.isArray(this.currentProtein.ligands) && this.currentProtein.ligands.length > 0) {
                    this.currentProtein.ligands.forEach(resn => {
                        const id = `lig-${resn}`;
                        const wrapper = document.createElement('div');
                        wrapper.className = 'checkbox-inline';
                        wrapper.innerHTML = `
                            <label class="checkbox-container">
                                <input type="checkbox" id="${id}" data-resn="${resn}">
                                <span class="checkmark"></span>
                                ${resn}
                            </label>`;
                        ligandContainer.appendChild(wrapper);
                    });
                } else {
                    ligandContainer.innerHTML = '<small>No ligands detected</small>';
                }
            }
        }
    }

    getSelectedChains() {
        const container = document.getElementById('chain-selection');
        if (!container) return [];
        return Array.from(container.querySelectorAll('input[type="checkbox"]:checked')).map(cb => cb.getAttribute('data-chain'));
    }

    getSelectedLigands() {
        const container = document.getElementById('ligand-selection');
        if (!container) return [];
        return Array.from(container.querySelectorAll('input[type="checkbox"]:checked')).map(cb => ({
            resn: cb.getAttribute('data-resn') || '',
            chain: cb.getAttribute('data-chain') || ''
        }));
    }

    previewPreparedStructure() {
        if (!this.preparedProtein) {
            alert('Please prepare a protein structure first');
            return;
        }

        // Show prepared structure preview
        document.getElementById('prepared-structure-preview').style.display = 'block';
        
        // Format removed components
        const removedText = this.preparedProtein.removed_components ? 
            Object.entries(this.preparedProtein.removed_components)
                .filter(([key, value]) => value > 0)
                .map(([key, value]) => `${key}: ${value}`)
                .join(', ') || 'None' : 'None';
        
        // Format added capping
        const addedText = this.preparedProtein.added_capping ? 
            Object.entries(this.preparedProtein.added_capping)
                .filter(([key, value]) => value > 0)
                .map(([key, value]) => `${key}: ${value}`)
                .join(', ') || 'None' : 'None';
        
        // Update structure info
        document.getElementById('original-atoms').textContent = this.preparedProtein.original_atoms.toLocaleString();
        document.getElementById('prepared-atoms').textContent = this.preparedProtein.prepared_atoms.toLocaleString();
        document.getElementById('removed-components').textContent = removedText;
        document.getElementById('added-capping').textContent = addedText;
        document.getElementById('preserved-ligands').textContent = this.preparedProtein.preserved_ligands;

        // Load 3D visualization of prepared structure
        this.loadPrepared3DVisualization();
    }

    downloadPreparedStructure() {
        if (!this.preparedProtein) {
            alert('Please prepare a structure first');
            return;
        }

        // Download prepared structure
        const blob = new Blob([this.preparedProtein.content], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `tleap_ready.pdb`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    downloadLigandFile() {
        if (!this.preparedProtein || !this.preparedProtein.ligand_present || !this.preparedProtein.ligand_content) {
            alert('No ligand file available. Please prepare structure with separate ligands enabled.');
            return;
        }

        // Download ligand file
        const ligandBlob = new Blob([this.preparedProtein.ligand_content], { type: 'text/plain' });
        const ligandUrl = URL.createObjectURL(ligandBlob);
        const ligandA = document.createElement('a');
        ligandA.href = ligandUrl;
        ligandA.download = `4_ligands_corrected.pdb`;
        document.body.appendChild(ligandA);
        ligandA.click();
        document.body.removeChild(ligandA);
        URL.revokeObjectURL(ligandUrl);
    }

    // 3D Visualization for prepared structure
    async loadPrepared3DVisualization() {
        if (!this.preparedProtein) return;

        try {
            // Initialize NGL stage for prepared structure if not already done
            if (!this.preparedNglStage) {
                this.preparedNglStage = new NGL.Stage("prepared-ngl-viewer", {
                    backgroundColor: "white",
                    quality: "medium"
                });
            }

            // Clear existing components
            this.preparedNglStage.removeAllComponents();

            // Create a blob from prepared PDB content
            const blob = new Blob([this.preparedProtein.content], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);

            // Load the prepared structure
            const component = await this.preparedNglStage.loadFile(url, {
                ext: "pdb",
                defaultRepresentation: false
            });

            // Add cartoon representation for each chain with consistent colors
            // This ensures each chain gets the same color as in Step 1
            // Use chains from parsed protein data (more reliable than structure API)
            const structureChains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                // Add representation for each chain that exists in the structure
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        component.addRepresentation("cartoon", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.9
                        });
                    }
                });
            } else {
                // Fallback: use chainid if color map not available
                component.addRepresentation("cartoon", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.9
                });
            }
            
            // Apply consistent chain colors after representation is added (backup)
            setTimeout(() => {
                this.applyConsistentChainColors(component);
            }, 500);

            // Add ball and stick for ligands (if any) - check for HETATM records
            component.addRepresentation("ball+stick", {
                sele: "hetero",
                color: "element",
                radius: 0.2,
                opacity: 0.8
            });

            // Auto-fit the view
            this.preparedNglStage.autoView();

            // Show controls
            document.getElementById('prepared-viewer-controls').style.display = 'flex';

            // Clean up the blob URL
            URL.revokeObjectURL(url);

        } catch (error) {
            console.error('Error loading prepared 3D visualization:', error);
        }
    }

    resetPreparedView() {
        if (this.preparedNglStage) {
            this.preparedNglStage.autoView();
        }
    }

    togglePreparedRepresentation() {
        if (!this.preparedNglStage) return;

        const components = this.preparedNglStage.compList;
        if (components.length === 0) return;

        const component = components[0];
        component.removeAllRepresentations();

        if (this.preparedRepresentation === 'cartoon') {
            // Switch to ball and stick
            component.addRepresentation("ball+stick", {
                color: "element",
                radius: 0.15
            });
            this.preparedRepresentation = 'ball+stick';
            document.getElementById('prepared-style-text').textContent = 'Ball & Stick';
        } else if (this.preparedRepresentation === 'ball+stick') {
            // Switch to surface with consistent chain colors
            // Use chains from parsed protein data (more reliable than structure API)
            const structureChains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        component.addRepresentation("surface", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.7
                        });
                    }
                });
            } else {
                component.addRepresentation("surface", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.7
                });
            }
            this.preparedRepresentation = 'surface';
            document.getElementById('prepared-style-text').textContent = 'Surface';
        } else {
            // Switch back to cartoon
            const chainColorFunc = this.getChainColorScheme(component);
            component.addRepresentation("cartoon", {
                sele: "protein",
                colorScheme: chainColorFunc,
                opacity: 0.8
            });

            // Add ball and stick for ligands
            if (this.preparedProtein.preserved_ligands !== 'None') {
                component.addRepresentation("ball+stick", {
                    sele: "hetero",
                    color: "element",
                    radius: 0.15
                });
            }

            this.preparedRepresentation = 'cartoon';
            document.getElementById('prepared-style-text').textContent = 'Mixed View';
        }
    }

    togglePreparedSpin() {
        if (!this.preparedNglStage) return;

        this.preparedIsSpinning = !this.preparedIsSpinning;
        this.preparedNglStage.setSpin(this.preparedIsSpinning);
    }

    // Missing Residues Methods
    async detectMissingResidues() {
        if (!this.currentProtein) {
            this.showMissingStatus('error', 'Please load a protein structure first');
            return;
        }

        const statusDiv = document.getElementById('missing-status');
        statusDiv.className = 'status-message info';
        statusDiv.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Detecting missing residues...';
        statusDiv.style.display = 'block';

        try {
            const response = await fetch('/api/detect-missing-residues', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                }
            });

            const result = await response.json();

            if (result.success) {
                this.missingResiduesInfo = result.missing_residues;
                this.missingResiduesPdbId = result.pdb_id;
                this.chainSequences = result.chain_sequences || {};

                if (Object.keys(result.missing_residues).length === 0) {
                    this.showMissingStatus('success', 'No missing residues detected in this structure!');
                    document.getElementById('missing-chains-section').style.display = 'none';
                    document.getElementById('trim-residues-section').style.display = 'none';
                    document.getElementById('build-complete-structure').disabled = true;
                    document.getElementById('missing-summary').style.display = 'none';
                } else {
                    this.showMissingStatus('success', `Found missing residues in ${Object.keys(result.missing_residues).length} chain(s)`);
                    this.renderMissingChains(result.chains_with_missing, result.missing_residues);
                    document.getElementById('missing-chains-section').style.display = 'block';
                    document.getElementById('missing-summary').style.display = 'block';
                    this.displayMissingSummary(result.missing_residues);
                    this.renderTrimControls(result.chains_with_missing);
                }
            } else {
                throw new Error(result.error || 'Failed to detect missing residues');
            }
        } catch (error) {
            console.error('Error detecting missing residues:', error);
            this.showMissingStatus('error', `Error: ${error.message}`);
        }
    }

    renderMissingChains(chainsWithMissing, missingResidues) {
        const container = document.getElementById('missing-chains-list');
        container.innerHTML = '';

        chainsWithMissing.forEach(chain => {
            const missingCount = missingResidues[chain]?.count || 0;
            const wrapper = document.createElement('div');
            wrapper.className = 'checkbox-inline';
            wrapper.innerHTML = `
                <label class="checkbox-container">
                    <input type="checkbox" id="missing-chain-${chain}" data-chain="${chain}" checked>
                    <span class="checkmark"></span>
                    Chain ${chain} (${missingCount} missing residues)
                </label>
            `;
            container.appendChild(wrapper);
        });

        // Update button state based on selections
        this.updateBuildButtonState();
        
        // Add event listeners to checkboxes
        container.querySelectorAll('input[type="checkbox"]').forEach(checkbox => {
            checkbox.addEventListener('change', () => this.updateBuildButtonState());
        });
    }

    updateBuildButtonState() {
        const container = document.getElementById('missing-chains-list');
        const selectedChains = Array.from(container.querySelectorAll('input[type="checkbox"]:checked'));
        const buildBtn = document.getElementById('build-complete-structure');
        const previewBtn = document.getElementById('preview-completed-structure');
        const downloadBtn = document.getElementById('download-completed-structure');

        if (selectedChains.length > 0) {
            buildBtn.disabled = false;
        } else {
            buildBtn.disabled = true;
            previewBtn.disabled = true;
            document.getElementById('preview-superimposed-structure').disabled = true;
            downloadBtn.disabled = true;
        }
    }

    displayMissingSummary(missingResidues) {
        const summaryContent = document.getElementById('missing-summary-content');
        let html = '';

        Object.entries(missingResidues).forEach(([chain, info]) => {
            html += `<div class="chain-missing-info">`;
            html += `<h4>Chain ${chain}: ${info.count} missing residues</h4>`;
            if (info.residues && info.residues.length > 0) {
                // Display all residues horizontally directly on green background
                const residueStrings = info.residues.map(
                    ([resname, resnum]) => `${resname} ${resnum}`
                );
                html += '<p class="missing-residues-horizontal">';
                html += residueStrings.join(', ');
                html += '</p>';
            }
            html += `</div>`;
        });

        summaryContent.innerHTML = html;
    }

    analyzeEdgeResidues(chain, missingResidues) {
        /**
         * Analyze missing residues to determine which ones are at edges
         * Returns: { n_terminal: {start, end, count}, c_terminal: {start, end, count} }
         */
        if (!missingResidues || !missingResidues[chain] || !missingResidues[chain].residues) {
            return { n_terminal: null, c_terminal: null };
        }

        const residues = missingResidues[chain].residues;
        if (residues.length === 0) {
            return { n_terminal: null, c_terminal: null };
        }

        // Extract residue numbers and sort
        const resNums = residues.map(([resname, resnum]) => resnum).sort((a, b) => a - b);
        
        // Find all consecutive ranges
        const ranges = [];
        if (resNums.length > 0) {
            let rangeStart = resNums[0];
            let rangeEnd = resNums[0];
            
            for (let i = 1; i < resNums.length; i++) {
                if (resNums[i] === rangeEnd + 1) {
                    // Consecutive, extend range
                    rangeEnd = resNums[i];
                } else {
                    // Gap found, save current range and start new one
                    ranges.push({ start: rangeStart, end: rangeEnd });
                    rangeStart = resNums[i];
                    rangeEnd = resNums[i];
                }
            }
            // Don't forget the last range
            ranges.push({ start: rangeStart, end: rangeEnd });
        }

        // Identify N-terminal edge (first range if it starts from a low number, typically 0 or 1)
        let nTerminal = null;
        if (ranges.length > 0) {
            const firstRange = ranges[0];
            // Consider it N-terminal if it starts from 0 or 1 (common in PDB files)
            if (firstRange.start <= 1) {
                nTerminal = {
                    start: firstRange.start,
                    end: firstRange.end,
                    count: firstRange.end - firstRange.start + 1
                };
            }
        }

        // Identify C-terminal edge (last range)
        // We consider it C-terminal if:
        // 1. There are multiple ranges (the last one is likely C-terminal since there are internal gaps)
        // 2. OR if there's only one range but it's not N-terminal (could be C-terminal only)
        let cTerminal = null;
        if (ranges.length > 0) {
            const lastRange = ranges[ranges.length - 1];
            
            if (ranges.length > 1) {
                // Multiple ranges - the last one is likely C-terminal (there are internal missing residues)
                cTerminal = {
                    start: lastRange.start,
                    end: lastRange.end,
                    count: lastRange.end - lastRange.start + 1
                };
            } else if (nTerminal === null) {
                // Only one range and it's not N-terminal
                // Check if sequence length is available to verify it's at the end
                const seqLength = this.chainSequences[chain] ? this.chainSequences[chain].length : 0;
                if (seqLength > 0 && lastRange.end >= seqLength - 5) {
                    // If the range ends within 5 residues of the sequence end, consider it C-terminal
                    cTerminal = {
                        start: lastRange.start,
                        end: lastRange.end,
                        count: lastRange.end - lastRange.start + 1
                    };
                } else if (lastRange.start > 10) {
                    // If the range starts well after the beginning (not N-terminal), 
                    // and there's only one range, it might be C-terminal
                    // But be conservative - only if it's clearly not at the start
                    cTerminal = {
                        start: lastRange.start,
                        end: lastRange.end,
                        count: lastRange.end - lastRange.start + 1
                    };
                }
            }
        }

        return {
            n_terminal: nTerminal,
            c_terminal: cTerminal
        };
    }

    updateTrimInfoBox(chainsWithMissing) {
        const infoBox = document.getElementById('trim-info-box-content');
        if (!infoBox || !this.missingResiduesInfo) {
            return;
        }

        let html = '<i class="fas fa-info-circle"></i> <strong>Note:</strong> ';
        
        const edgeInfo = [];
        chainsWithMissing.forEach(chain => {
            const edges = this.analyzeEdgeResidues(chain, this.missingResiduesInfo);
            const chainInfo = [];
            
            if (edges.n_terminal) {
                chainInfo.push(`residues ${edges.n_terminal.start}-${edges.n_terminal.end} from N-terminal`);
            }
            if (edges.c_terminal) {
                chainInfo.push(`residues ${edges.c_terminal.start}-${edges.c_terminal.end} from C-terminal`);
            }
            
            if (chainInfo.length > 0) {
                edgeInfo.push(`Chain ${chain}: ${chainInfo.join(' and ')}`);
            }
        });

        if (edgeInfo.length > 0) {
            html += 'Only missing residues at the edges can be trimmed. ';
            html += edgeInfo.join('; ') + '. ';
            html += 'Missing residues in internal loops (discontinuities in the middle) cannot be trimmed and will be filled by ESMFold.';
        } else {
            html += 'Only missing residues at the N-terminal edge (beginning) and C-terminal edge (end) can be trimmed. ';
            html += 'Missing residues in internal loops (discontinuities in the middle of the sequence) cannot be trimmed using this tool and will be filled by ESMFold.';
        }

        infoBox.innerHTML = html;
    }

    renderTrimControls(chainsWithMissing) {
        const container = document.getElementById('trim-residues-list');
        container.innerHTML = '';

        if (!this.chainSequences || Object.keys(this.chainSequences).length === 0) {
            return;
        }

        // Update the info box with dynamic information
        this.updateTrimInfoBox(chainsWithMissing);

        chainsWithMissing.forEach(chain => {
            const sequence = this.chainSequences[chain] || '';
            const seqLength = sequence.length;
            
            // Get edge residue information for this chain
            const edges = this.analyzeEdgeResidues(chain, this.missingResiduesInfo);
            
            // Calculate max values based on detected edge residues
            const nTerminalMax = edges.n_terminal ? edges.n_terminal.count : 0;
            const cTerminalMax = edges.c_terminal ? edges.c_terminal.count : 0;
            
            // Build N-terminal label with limit info
            let nTerminalLabel = 'N-terminal:';
            if (edges.n_terminal) {
                nTerminalLabel += ` <span class="trim-limit">(max: ${nTerminalMax})</span>`;
            } else {
                nTerminalLabel += ` <span class="trim-limit" style="color: #6c757d; font-style: italic;">(no missing residues)</span>`;
            }
            
            // Build C-terminal label with limit info
            let cTerminalLabel = 'C-terminal:';
            if (edges.c_terminal) {
                cTerminalLabel += ` <span class="trim-limit">(max: ${cTerminalMax})</span>`;
            } else {
                cTerminalLabel += ` <span class="trim-limit" style="color: #6c757d; font-style: italic;">(no missingresidues)</span>`;
            }

            const wrapper = document.createElement('div');
            wrapper.className = 'trim-chain-controls';
            wrapper.innerHTML = `
                <h5>Chain ${chain} (${seqLength} residues)</h5>
                <div class="trim-inputs">
                    <div class="trim-input-group">
                        <label>${nTerminalLabel}</label>
                        <input type="number" 
                               id="trim-n-${chain}" 
                               data-chain="${chain}" 
                               min="0" 
                               max="${nTerminalMax}" 
                               value="0"
                               class="trim-n-input"
                               ${nTerminalMax > 0 ? `data-max-edge="${nTerminalMax}"` : 'disabled'}
                               ${nTerminalMax === 0 ? 'style="background-color: #e9ecef; cursor: not-allowed;"' : ''}>
                        <span>residues</span>
                    </div>
                    <div class="trim-input-group">
                        <label>${cTerminalLabel}</label>
                        <input type="number" 
                               id="trim-c-${chain}" 
                               data-chain="${chain}" 
                               min="0" 
                               max="${cTerminalMax}" 
                               value="0"
                               class="trim-c-input"
                               ${cTerminalMax > 0 ? `data-max-edge="${cTerminalMax}"` : 'disabled'}
                               ${cTerminalMax === 0 ? 'style="background-color: #e9ecef; cursor: not-allowed;"' : ''}>
                        <span>residues</span>
                    </div>
                </div>
                <div class="trim-info" id="trim-info-${chain}">
                    Original length: ${seqLength} residues
                </div>
            `;
            container.appendChild(wrapper);

            // Add event listeners to update info and enforce limits
            const nInput = wrapper.querySelector(`#trim-n-${chain}`);
            const cInput = wrapper.querySelector(`#trim-c-${chain}`);
            const infoDiv = wrapper.querySelector(`#trim-info-${chain}`);

            // Enforce max limits based on edge residues (only if there are edge residues)
            if (nTerminalMax > 0) {
                nInput.addEventListener('input', () => {
                    const value = parseInt(nInput.value) || 0;
                    if (value > nTerminalMax) {
                        nInput.value = nTerminalMax;
                    }
                });
            } else {
                // Disable input if no edge residues
                nInput.disabled = true;
            }
            
            if (cTerminalMax > 0) {
                cInput.addEventListener('input', () => {
                    const value = parseInt(cInput.value) || 0;
                    if (value > cTerminalMax) {
                        cInput.value = cTerminalMax;
                    }
                });
            } else {
                // Disable input if no edge residues
                cInput.disabled = true;
            }

            const updateInfo = () => {
                const nTrim = parseInt(nInput.value) || 0;
                const cTrim = parseInt(cInput.value) || 0;
                const totalTrim = nTrim + cTrim;
                const newLength = seqLength - totalTrim;
                
                // Check if values exceed edge limits
                let warningMsg = '';
                if (nTerminalMax > 0 && nTrim > nTerminalMax) {
                    warningMsg = `<span style="color: #dc3545;">Warning: N-terminal trim (${nTrim}) exceeds edge limit (${nTerminalMax})</span>`;
                } else if (cTerminalMax > 0 && cTrim > cTerminalMax) {
                    warningMsg = `<span style="color: #dc3545;">Warning: C-terminal trim (${cTrim}) exceeds edge limit (${cTerminalMax})</span>`;
                } else if (nTerminalMax === 0 && nTrim > 0) {
                    warningMsg = `<span style="color: #dc3545;">Warning: No N-terminal edge residues to trim</span>`;
                } else if (cTerminalMax === 0 && cTrim > 0) {
                    warningMsg = `<span style="color: #dc3545;">Warning: No C-terminal edge residues to trim</span>`;
                } else if (totalTrim >= seqLength) {
                    warningMsg = `<span style="color: #dc3545;">Error: Total trim (${totalTrim}) exceeds sequence length (${seqLength})</span>`;
                } else if (newLength <= 0) {
                    warningMsg = `<span style="color: #dc3545;">Error: Resulting sequence would be empty</span>`;
                } else {
                    let infoText = `Original: ${seqLength} residues → Trimmed: ${newLength} residues (removing ${nTrim} from N-term, ${cTrim} from C-term)`;
                    if (nTerminalMax > 0 || cTerminalMax > 0) {
                        infoText += `<br><small style="color: #6c757d;">Edge limits: N-term max ${nTerminalMax}, C-term max ${cTerminalMax}</small>`;
                    } else {
                        infoText += `<br><small style="color: #6c757d;">No edge residues available for trimming</small>`;
                    }
                    infoDiv.innerHTML = infoText;
                }
                
                if (warningMsg) {
                    infoDiv.innerHTML = warningMsg;
                }
            };

            nInput.addEventListener('input', updateInfo);
            cInput.addEventListener('input', updateInfo);
            updateInfo(); // Initial update
        });

        // Show the trim section
        document.getElementById('trim-residues-section').style.display = 'block';
    }

    async applyTrimming() {
        if (!this.chainSequences || !this.missingResiduesPdbId) {
            this.showTrimStatus('error', 'No chain sequences available. Please detect missing residues first.');
            return;
        }

        const trimStatusDiv = document.getElementById('trim-status');
        trimStatusDiv.className = 'status-message info';
        trimStatusDiv.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Applying trimming...';
        trimStatusDiv.style.display = 'block';

        try {
            // Collect trim specifications from inputs
            const trimSpecs = {};
            const nInputs = document.querySelectorAll('.trim-n-input');
            const cInputs = document.querySelectorAll('.trim-c-input');

            nInputs.forEach(nInput => {
                const chain = nInput.getAttribute('data-chain');
                const nTrim = parseInt(nInput.value) || 0;
                const cInput = document.querySelector(`#trim-c-${chain}`);
                const cTrim = parseInt(cInput.value) || 0;

                if (nTrim > 0 || cTrim > 0) {
                    trimSpecs[chain] = {
                        n_terminal: nTrim,
                        c_terminal: cTrim
                    };
                }
            });

            if (Object.keys(trimSpecs).length === 0) {
                this.showTrimStatus('info', 'No trimming specified. Enter values > 0 to trim residues.');
                return;
            }

            // Call API to apply trimming
            const response = await fetch('/api/trim-residues', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    pdb_id: this.missingResiduesPdbId,
                    chain_sequences: this.chainSequences,
                    trim_specs: trimSpecs
                })
            });

            const result = await response.json();

            if (result.success) {
                // Update stored sequences with trimmed versions
                this.chainSequences = result.trimmed_sequences;
                
                // Update trim info displays
                Object.entries(result.trim_info).forEach(([chain, info]) => {
                    const infoDiv = document.getElementById(`trim-info-${chain}`);
                    if (infoDiv) {
                        infoDiv.innerHTML = `
                            <strong>Trimmed!</strong> Original: ${info.original_length} → 
                            Trimmed: ${info.trimmed_length} residues 
                            (removed ${info.n_terminal_trimmed} from N-term, ${info.c_terminal_trimmed} from C-term)
                        `;
                        infoDiv.style.color = '#155724';
                    }
                });

                this.showTrimStatus('success', result.message);
            } else {
                throw new Error(result.error || 'Failed to apply trimming');
            }
        } catch (error) {
            console.error('Error applying trimming:', error);
            this.showTrimStatus('error', `Error: ${error.message}`);
        }
    }

    showTrimStatus(type, message) {
        const statusDiv = document.getElementById('trim-status');
        statusDiv.className = `status-message ${type}`;
        statusDiv.textContent = message;
        statusDiv.style.display = 'block';

        if (type === 'success') {
            setTimeout(() => {
                statusDiv.style.display = 'none';
            }, 5000);
        }
    }

    async buildCompletedStructure() {
        const container = document.getElementById('missing-chains-list');
        const selectedChains = Array.from(container.querySelectorAll('input[type="checkbox"]:checked'))
            .map(cb => cb.getAttribute('data-chain'));

        if (selectedChains.length === 0) {
            this.showMissingStatus('error', 'Please select at least one chain to complete');
            return;
        }

        const buildBtn = document.getElementById('build-complete-structure');
        const originalText = buildBtn.innerHTML;
        buildBtn.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Building...';
        buildBtn.disabled = true;

        try {
            // Prepare request body with optional trimmed sequences
            const requestBody = {
                selected_chains: selectedChains
            };
            
            // Include trimmed sequences if available (they may have been trimmed)
            if (this.chainSequences && Object.keys(this.chainSequences).length > 0) {
                // Only include sequences for selected chains
                const selectedSequences = {};
                selectedChains.forEach(chain => {
                    if (this.chainSequences[chain]) {
                        selectedSequences[chain] = this.chainSequences[chain];
                    }
                });
                if (Object.keys(selectedSequences).length > 0) {
                    requestBody.chain_sequences = selectedSequences;
                }
            }

            const response = await fetch('/api/build-completed-structure', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(requestBody)
            });

            const result = await response.json();

            if (result.success) {
                this.completedProtein = {
                    content: result.completed_structure,
                    completed_chains: result.completed_chains
                };

                this.showMissingStatus('success', result.message);
                document.getElementById('preview-completed-structure').disabled = false;
                document.getElementById('preview-superimposed-structure').disabled = false;
                document.getElementById('download-completed-structure').disabled = false;
            } else {
                throw new Error(result.error || 'Failed to build completed structure');
            }
        } catch (error) {
            console.error('Error building completed structure:', error);
            this.showMissingStatus('error', `Error: ${error.message}`);
        } finally {
            buildBtn.innerHTML = originalText;
            buildBtn.disabled = false;
        }
    }

    async previewCompletedStructure() {
        if (!this.completedProtein) {
            // Try to fetch from server
            try {
                const response = await fetch('/api/get-completed-structure');
                const result = await response.json();
                
                if (result.success && result.exists) {
                    this.completedProtein = {
                        content: result.content
                    };
                } else {
                    alert('Completed structure not found. Please build it first.');
                    return;
                }
            } catch (error) {
                alert('Error loading completed structure: ' + error.message);
                return;
            }
        }

        // Get original structure
        if (!this.currentProtein) {
            alert('Original structure not found. Please load a PDB file first.');
            return;
        }

        // Show preview in the same tab
        const previewDiv = document.getElementById('completed-structure-preview');
        previewDiv.style.display = 'block';

        // Load both structures side by side
        try {
            // Load original structure
            await this.loadOriginalStructureViewer();
            
            // Load completed structure
            await this.loadCompletedStructureViewer();
        } catch (error) {
            console.error('Error previewing structures:', error);
            alert('Error loading 3D visualization: ' + error.message);
        }
    }

    async loadOriginalStructureViewer() {
        try {
            // Initialize NGL stage for original structure if not already done
            if (!this.originalNglStage) {
                this.originalNglStage = new NGL.Stage("original-ngl-viewer", {
                    backgroundColor: "white",
                    quality: "medium"
                });
            }

            // Clear existing components
            this.originalNglStage.removeAllComponents();

            // Create a blob from original PDB content
            const blob = new Blob([this.currentProtein.content], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);

            // Load the original structure
            const component = await this.originalNglStage.loadFile(url, {
                ext: "pdb",
                defaultRepresentation: false
            });

            // Add cartoon representation for each chain with consistent colors
            // This ensures each chain gets the same color as in Step 1
            // Use chains from parsed protein data (more reliable than structure API)
            const structureChains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                // Add representation for each chain that exists in the structure
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        component.addRepresentation("cartoon", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.9
                        });
                    }
                });
            } else {
                // Fallback: use chainid if color map not available
                component.addRepresentation("cartoon", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.9
                });
            }
            
            // Apply consistent chain colors after representation is added (backup)
            setTimeout(() => {
                this.applyConsistentChainColors(component);
            }, 500);

            // Add ball and stick for ligands if present
            if (this.currentProtein.ligands && this.currentProtein.ligands.length > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "hetero",
                    color: "element",
                    radius: 0.15
                });
            }

            // Auto-fit the view
            this.originalNglStage.autoView();

            // Show controls
            document.getElementById('original-viewer-controls').style.display = 'flex';

            // Clean up the blob URL
            URL.revokeObjectURL(url);
        } catch (error) {
            console.error('Error loading original structure viewer:', error);
            throw error;
        }
    }

    async previewSuperimposedStructure() {
        // Show preview section
        const previewDiv = document.getElementById('superimposed-structure-preview');
        previewDiv.style.display = 'block';

        try {
            // Try to use in-memory data first (faster)
            let originalText, completedText;
            
            if (this.currentProtein && this.currentProtein.content) {
                // Use already loaded original structure
                originalText = this.currentProtein.content;
            } else {
                // Fallback: fetch from server
                const originalResponse = await fetch('/api/get-file?filename=0_original_input.pdb');
                if (!originalResponse.ok) {
                    throw new Error('Failed to load original structure file.');
                }
                originalText = await originalResponse.text();
            }
            
            if (this.completedProtein && this.completedProtein.content) {
                // Use already loaded completed structure
                completedText = this.completedProtein.content;
            } else {
                // Fallback: fetch from server
                const completedResponse = await fetch('/api/get-file?filename=0_complete_structure.pdb');
                if (!completedResponse.ok) {
                    throw new Error('Failed to load completed structure file.');
                }
                completedText = await completedResponse.text();
            }

            // Initialize NGL stage for superimposed view if not already done
            if (!this.superimposedNglStage) {
                this.superimposedNglStage = new NGL.Stage("superimposed-ngl-viewer", {
                    backgroundColor: "white",
                    quality: "medium"
                });
            }

            // Clear existing components
            this.superimposedNglStage.removeAllComponents();

            // Create blobs from PDB content
            const originalBlob = new Blob([originalText], { type: 'text/plain' });
            const completedBlob = new Blob([completedText], { type: 'text/plain' });
            const originalUrl = URL.createObjectURL(originalBlob);
            const completedUrl = URL.createObjectURL(completedBlob);

            // Load original structure with original colors
            const originalComponent = await this.superimposedNglStage.loadFile(originalUrl, {
                ext: "pdb",
                defaultRepresentation: false
            });

            // Apply original chain colors if available
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                const structureChains = this.currentProtein && this.currentProtein.chains ? this.currentProtein.chains : [];
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        originalComponent.addRepresentation("cartoon", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.8
                        });
                    }
                });
            } else {
                // Fallback: use chainid
                originalComponent.addRepresentation("cartoon", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.8
                });
            }

            // Load completed structure with different colors
            const completedComponent = await this.superimposedNglStage.loadFile(completedUrl, {
                ext: "pdb",
                defaultRepresentation: false
            });

            // Get chains from completed structure - parse from PDB content or use original chains
            const completedChains = [];
            
            // Method 1: Parse chain IDs from PDB content
            const chainSet = new Set();
            const lines = completedText.split('\n');
            for (const line of lines) {
                if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
                    const chainId = line.charAt(21); // Chain ID is at position 21
                    if (chainId && chainId.trim() !== '') {
                        chainSet.add(chainId);
                    }
                }
            }
            completedChains.push(...Array.from(chainSet).sort());
            
            // Method 2: Fallback to original structure chains if parsing didn't work
            if (completedChains.length === 0 && this.currentProtein && this.currentProtein.chains) {
                completedChains.push(...this.currentProtein.chains);
            }
            
            // Use a different color scheme for completed structure
            const completedColors = {
                'A': 'orange',
                'B': 'purple',
                'C': 'pink',
                'D': 'cyan',
                'E': 'yellow',
                'F': 'magenta',
                'G': 'lime',
                'H': 'coral'
            };

            if (completedChains.length > 0) {
                completedChains.forEach((chain, index) => {
                    // Use different colors for completed structure
                    const color = completedColors[chain] || `hsl(${(index * 60) % 360}, 70%, 50%)`;
                    completedComponent.addRepresentation("cartoon", {
                        sele: `:${chain}`,
                        color: color,
                        opacity: 0.7
                    });
                });
            } else {
                // Fallback: use chainid color scheme if we can't determine chains
                completedComponent.addRepresentation("cartoon", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.7
                });
            }

            // Store completed chains for later use in toggle functions
            this.completedChains = completedChains;

            // Auto-fit the view
            this.superimposedNglStage.autoView();

            // Show controls
            document.getElementById('superimposed-viewer-controls').style.display = 'flex';

            // Clean up blob URLs
            URL.revokeObjectURL(originalUrl);
            URL.revokeObjectURL(completedUrl);

            // Store components for control functions
            this.superimposedOriginalComponent = originalComponent;
            this.superimposedCompletedComponent = completedComponent;
            this.superimposedRepresentationType = 'cartoon';
            this.superimposedIsSpinning = false;

        } catch (error) {
            console.error('Error loading superimposed structures:', error);
            alert('Error loading superimposed visualization: ' + error.message);
        }
    }

    resetSuperimposedView() {
        if (this.superimposedNglStage) {
            this.superimposedNglStage.autoView();
        }
    }

    toggleSuperimposedRepresentation() {
        if (!this.superimposedOriginalComponent || !this.superimposedCompletedComponent) {
            return;
        }

        // Remove existing representations
        this.superimposedOriginalComponent.removeAllRepresentations();
        this.superimposedCompletedComponent.removeAllRepresentations();

        const styleText = document.getElementById('superimposed-style-text');
        
        if (this.superimposedRepresentationType === 'cartoon') {
            // Switch to surface
            this.superimposedRepresentationType = 'surface';
            styleText.textContent = 'Surface';
            
            // Original structure
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                const structureChains = this.currentProtein && this.currentProtein.chains ? this.currentProtein.chains : [];
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        this.superimposedOriginalComponent.addRepresentation("surface", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.6
                        });
                    }
                });
            } else {
                this.superimposedOriginalComponent.addRepresentation("surface", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.6
                });
            }

            // Completed structure - use stored chains or original chains
            let completedChains = [];
            
            // Try to get from stored completed chains (if we stored them)
            if (this.completedChains && this.completedChains.length > 0) {
                completedChains = this.completedChains;
            } else if (this.currentProtein && this.currentProtein.chains) {
                // Fallback to original chains
                completedChains = this.currentProtein.chains;
            }
            
            const completedColors = {
                'A': 'orange',
                'B': 'purple',
                'C': 'pink',
                'D': 'cyan',
                'E': 'yellow',
                'F': 'magenta',
                'G': 'lime',
                'H': 'coral'
            };
            
            if (completedChains.length > 0) {
                completedChains.forEach((chain, index) => {
                    const color = completedColors[chain] || `hsl(${(index * 60) % 360}, 70%, 50%)`;
                    this.superimposedCompletedComponent.addRepresentation("surface", {
                        sele: `:${chain}`,
                        color: color,
                        opacity: 0.5
                    });
                });
            } else {
                // Fallback: use chainid color scheme
                this.superimposedCompletedComponent.addRepresentation("surface", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.5
                });
            }
        } else {
            // Switch to cartoon
            this.superimposedRepresentationType = 'cartoon';
            styleText.textContent = 'Cartoon';
            
            // Original structure
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                const structureChains = this.currentProtein && this.currentProtein.chains ? this.currentProtein.chains : [];
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        this.superimposedOriginalComponent.addRepresentation("cartoon", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.8
                        });
                    }
                });
            } else {
                this.superimposedOriginalComponent.addRepresentation("cartoon", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.8
                });
            }

            // Completed structure - use stored chains or original chains
            let completedChains = [];
            
            // Try to get from stored completed chains (if we stored them)
            if (this.completedChains && this.completedChains.length > 0) {
                completedChains = this.completedChains;
            } else if (this.currentProtein && this.currentProtein.chains) {
                // Fallback to original chains
                completedChains = this.currentProtein.chains;
            }
            
            const completedColors = {
                'A': 'orange',
                'B': 'purple',
                'C': 'pink',
                'D': 'cyan',
                'E': 'yellow',
                'F': 'magenta',
                'G': 'lime',
                'H': 'coral'
            };
            
            if (completedChains.length > 0) {
                completedChains.forEach((chain, index) => {
                    const color = completedColors[chain] || `hsl(${(index * 60) % 360}, 70%, 50%)`;
                    this.superimposedCompletedComponent.addRepresentation("cartoon", {
                        sele: `:${chain}`,
                        color: color,
                        opacity: 0.7
                    });
                });
            } else {
                // Fallback: use chainid color scheme
                this.superimposedCompletedComponent.addRepresentation("cartoon", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.7
                });
            }
        }
    }

    toggleSuperimposedSpin() {
        if (!this.superimposedNglStage) {
            return;
        }
        this.superimposedIsSpinning = !this.superimposedIsSpinning;
        this.superimposedNglStage.setSpin(this.superimposedIsSpinning);
    }

    async loadCompletedStructureViewer() {
        try {
            // Initialize NGL stage for completed structure if not already done
            if (!this.completedNglStage) {
                this.completedNglStage = new NGL.Stage("completed-ngl-viewer", {
                    backgroundColor: "white",
                    quality: "medium"
                });
            }

            // Clear existing components
            this.completedNglStage.removeAllComponents();

            // Create a blob from completed PDB content
            const blob = new Blob([this.completedProtein.content], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);

            // Load the completed structure
            const component = await this.completedNglStage.loadFile(url, {
                ext: "pdb",
                defaultRepresentation: false
            });

            // Add cartoon representation for each chain with consistent colors
            // This ensures each chain gets the same color as in Step 1
            // Use chains from parsed protein data (more reliable than structure API)
            const structureChains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                // Add representation for each chain that exists in the structure
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        component.addRepresentation("cartoon", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.9
                        });
                    }
                });
            } else {
                // Fallback: use chainid if color map not available
                component.addRepresentation("cartoon", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.9
                });
            }
            
            // Apply consistent chain colors after representation is added (backup)
            setTimeout(() => {
                this.applyConsistentChainColors(component);
            }, 500);

            // Add ball and stick for ligands if present
            if (this.currentProtein.ligands && this.currentProtein.ligands.length > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "hetero",
                    color: "element",
                    radius: 0.15
                });
            }

            // Auto-fit the view
            this.completedNglStage.autoView();

            // Show controls
            document.getElementById('completed-viewer-controls').style.display = 'flex';

            // Clean up the blob URL
            URL.revokeObjectURL(url);
        } catch (error) {
            console.error('Error loading completed structure viewer:', error);
            throw error;
        }
    }

    resetCompletedView() {
        if (this.completedNglStage) {
            this.completedNglStage.autoView();
        }
    }

    toggleCompletedRepresentation() {
        if (!this.completedNglStage) return;

        const components = this.completedNglStage.compList;
        if (components.length === 0) return;

        const component = components[0];
        component.removeAllRepresentations();

        if (this.completedRepresentation === 'cartoon') {
            // Switch to ball and stick
            component.addRepresentation("ball+stick", {
                color: "element",
                radius: 0.15
            });
            this.completedRepresentation = 'ball+stick';
            document.getElementById('completed-style-text').textContent = 'Ball & Stick';
        } else if (this.completedRepresentation === 'ball+stick') {
            // Switch to surface with consistent chain colors
            // Use chains from parsed protein data (more reliable than structure API)
            const structureChains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        component.addRepresentation("surface", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.7
                        });
                    }
                });
            } else {
                component.addRepresentation("surface", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.7
                });
            }
            this.completedRepresentation = 'surface';
            document.getElementById('completed-style-text').textContent = 'Surface';
        } else {
            // Switch back to cartoon
            const chainColorFunc = this.getChainColorScheme(component);
            component.addRepresentation("cartoon", {
                sele: "protein",
                colorScheme: chainColorFunc,
                opacity: 0.8
            });
            // Add ligands if present
            if (this.currentProtein && this.currentProtein.ligands && this.currentProtein.ligands.length > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "hetero",
                    color: "element",
                    radius: 0.15
                });
            }
            this.completedRepresentation = 'cartoon';
            document.getElementById('completed-style-text').textContent = 'Mixed';
        }
    }

    toggleCompletedSpin() {
        if (!this.completedNglStage) return;

        this.completedIsSpinning = !this.completedIsSpinning;
        this.completedNglStage.setSpin(this.completedIsSpinning);
    }

    // Original structure viewer controls
    resetOriginalView() {
        if (this.originalNglStage) {
            this.originalNglStage.autoView();
        }
    }

    toggleOriginalRepresentation() {
        if (!this.originalNglStage) return;

        const components = this.originalNglStage.compList;
        if (components.length === 0) return;

        const component = components[0];
        component.removeAllRepresentations();

        if (this.originalRepresentation === 'cartoon') {
            // Switch to ball and stick
            component.addRepresentation("ball+stick", {
                color: "element",
                radius: 0.15
            });
            this.originalRepresentation = 'ball+stick';
            document.getElementById('original-style-text').textContent = 'Ball & Stick';
        } else if (this.originalRepresentation === 'ball+stick') {
            // Switch to surface with consistent chain colors
            // Use chains from parsed protein data (more reliable than structure API)
            const structureChains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            if (this.chainColorMap && Object.keys(this.chainColorMap).length > 0) {
                structureChains.forEach((chain) => {
                    if (this.chainColorMap[chain]) {
                        component.addRepresentation("surface", {
                            sele: `:${chain}`,
                            color: this.chainColorMap[chain],
                            opacity: 0.7
                        });
                    }
                });
            } else {
                component.addRepresentation("surface", {
                    sele: "protein",
                    colorScheme: "chainid",
                    opacity: 0.7
                });
            }
            this.originalRepresentation = 'surface';
            document.getElementById('original-style-text').textContent = 'Surface';
        } else {
            // Switch back to cartoon
            const chainColorFunc = this.getChainColorScheme(component);
            component.addRepresentation("cartoon", {
                sele: "protein",
                colorScheme: chainColorFunc,
                opacity: 0.8
            });
            if (this.currentProtein.ligands && this.currentProtein.ligands.length > 0) {
                component.addRepresentation("ball+stick", {
                    sele: "hetero",
                    color: "element",
                    radius: 0.15
                });
            }
            this.originalRepresentation = 'cartoon';
            document.getElementById('original-style-text').textContent = 'Mixed';
        }
    }

    toggleOriginalSpin() {
        if (!this.originalNglStage) return;

        this.originalIsSpinning = !this.originalIsSpinning;
        this.originalNglStage.setSpin(this.originalIsSpinning);
    }


    downloadCompletedStructure() {
        if (!this.completedProtein) {
            alert('Completed structure not found. Please build it first.');
            return;
        }

        const blob = new Blob([this.completedProtein.content], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = '0_complete_structure.pdb';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    showMissingStatus(type, message) {
        const statusDiv = document.getElementById('missing-status');
        statusDiv.className = `status-message ${type}`;
        statusDiv.textContent = message;
        statusDiv.style.display = 'block';

        if (type === 'success') {
            setTimeout(() => {
                statusDiv.style.display = 'none';
            }, 5000);
        }
    }

    applyConsistentChainColors(component) {
        // Apply consistent colors to chains using NGL's API
        if (!component || !component.structure || !this.chainColorMap || Object.keys(this.chainColorMap).length === 0) {
            console.warn('Cannot apply chain colors: missing component, structure, or color map');
            return;
        }
        
        try {
            // Get all chains - use parsed protein data (more reliable than structure API)
            const chains = (this.currentProtein && this.currentProtein.chains) ? this.currentProtein.chains : [];
            console.log('Applying colors to chains:', chains, 'Color map:', this.chainColorMap);
            
            // Apply colors to each representation using setColorByChain
            component.reprList.forEach((repr) => {
                if (repr.type === 'cartoon' || repr.type === 'surface') {
                    chains.forEach((chain) => {
                        if (this.chainColorMap[chain]) {
                            try {
                                const color = this.chainColorMap[chain];
                                // Use setColorByChain if available, otherwise use setColor
                                if (repr.setColorByChain) {
                                    repr.setColorByChain(color, chain);
                                } else {
                                    // Fallback: use setColor with chain selection
                                    repr.setColor(color, `chain ${chain}`);
                                }
                                console.log(`Applied color ${color} to chain ${chain}`);
                            } catch (err) {
                                console.warn(`Could not set color for chain ${chain}:`, err);
                            }
                        }
                    });
                }
            });
        } catch (error) {
            console.warn('Could not apply consistent chain colors:', error);
        }
    }
}

// Initialize the application when the page loads
function initializeApp() {
    console.log('Initializing mdPipeline...'); // Debug log
    window.mdPipeline = new MDSimulationPipeline();
    console.log('mdPipeline initialized:', window.mdPipeline); // Debug log
}

// Try to initialize immediately if DOM is already loaded
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initializeApp);
} else {
    // DOM is already loaded
    initializeApp();
}

// Add some utility functions for better UX
function formatNumber(num) {
    return num.toLocaleString();
}

function formatTime(seconds) {
    const hours = Math.floor(seconds / 3600);
    const minutes = Math.floor((seconds % 3600) / 60);
    const secs = seconds % 60;
    return `${hours}h ${minutes}m ${secs}s`;
}
