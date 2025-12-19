// MD Simulation Pipeline JavaScript
console.log('Script loading...'); // Debug log

class MDSimulationPipeline {
    constructor() {
        this.currentProtein = null;
        this.preparedProtein = null;
        this.simulationParams = {};
        this.generatedFiles = {};
        this.nglStage = null;
        this.preparedNglStage = null;
        this.currentRepresentation = 'cartoon';
        this.preparedRepresentation = 'cartoon';
        this.isSpinning = false;
        this.preparedIsSpinning = false;
        this.currentTabIndex = 0;
        this.tabOrder = ['protein-loading', 'structure-prep', 'simulation-params', 'simulation-steps', 'file-generation'];
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
        try {
            // Clean output folder when new PDB is loaded
            try {
                await fetch('/api/clean-output', { method: 'POST' });
            } catch (error) {
                console.log('Could not clean output folder:', error);
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

        document.getElementById('protein-preview').style.display = 'block';
        
        // Load 3D visualization
        this.load3DVisualization();

        // Also refresh chain/ligand lists when protein info is displayed
        this.renderChainAndLigandSelections();
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

            // Add cartoon representation for protein with chain-based colors
            component.addRepresentation("cartoon", {
                sele: "protein",
                colorScheme: "chainname",
                opacity: 0.9
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
            component.addRepresentation("surface", {
                sele: "protein",
                colorScheme: "chainname",
                opacity: 0.7
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

            // Add cartoon representation for protein with chain-based colors
            component.addRepresentation("cartoon", {
                sele: "protein",
                colorScheme: "chainname",
                opacity: 0.9
            });

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
            // Switch to surface
            component.addRepresentation("surface", {
                sele: "protein",
                colorScheme: "chainname",
                opacity: 0.7
            });
            this.preparedRepresentation = 'surface';
            document.getElementById('prepared-style-text').textContent = 'Surface';
        } else {
            // Switch back to cartoon
            component.addRepresentation("cartoon", {
                sele: "protein",
                colorScheme: "chainname",
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
