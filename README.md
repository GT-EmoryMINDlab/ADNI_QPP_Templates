<p align="center">
  <h1>Widespread Spatiotemporal Patterns of Functional Brain Networks in Longitudinal Progression of Alzheimer’s Disease</h1>
</p>

<p align="center">
  <strong>Theodore J. LaGrow<sup>1</sup>, Vaibhavi Itkyal<sup>2</sup>, Harrison Watters<sup>4</sup>, Kyle M. Jensen<sup>2</sup>, Ram Ballem<sup>2</sup>, Wen-ju Pan<sup>3</sup>, Armin Iraji<sup>2</sup>, Vince D. Calhoun<sup>1,2,3</sup>, Shella Keilholz<sup>3,4</sup></strong>
</p>

<p align="center">
  <sup>1</sup>School of Electrical and Computer Engineering, Georgia Institute of Technology <br>
  <sup>2</sup>Tri-institutional Center for Translational Research in Neuroimaging and Data Science, GSU, GT, and Emory University <br>
  <sup>3</sup>Wallace H. Coulter Department of Biomedical Engineering, Georgia Institute of Technology & Emory University <br>
  <sup>4</sup>Department of Neuroscience, Emory University
</p>

<p align="center">
  <img src="figures/ADNI_FULL_QPP_visualization.gif" alt="Visualization">

</p>

<p align="center">
  Repository corresponding to EMBC 2025 submission.
</p>

## Abstract
Alzheimer’s Disease (AD) is characterized by progressive functional network disruptions that precede cognitive decline, yet traditional functional connectivity analyses often fail to capture transient network instabilities critical for early diagnosis. This study investigates the role of Quasi-Periodic Patterns (QPPs) in identifying disease-related connectivity changes across longitudinal stable disease stages (sNC, sMCI, sDAT) and transitioning (uNC, pMCI) AD cohorts using resting-state fMRI data from the Alzheimer’s Disease Neuroimaging Initiative. By integrating QPP occurrences with intrinsic connectivity networks (ICNs), we assessed network integrity across disease stages, with statistical significance evaluated using the Kruskal-Wallis test and Dunn’s test for post-hoc analysis. Results revealed a progressive decline in functional connectivity integrity, with early impairments in subcortical and executive function networks in stable groups, followed by widespread disconnection in higher cognition, sensorimotor, and visual networks at later stages. Transitioning AD groups exhibited earlier disruptions in visual and cerebellar networks, suggesting their potential as early biomarkers for disease onset. The occurrence of QPPs decreased significantly with disease progression, indicating an increase in functional disconnection. These findings highlight the synergy between QPPs and ICNs as a dynamic and sensitive biomarker framework for AD progression. Future research should further explore this integration within multimodal imaging and clinical diagnostic frameworks to enhance early detection and intervention strategies.

## Data Access
The Alzheimer’s Disease Neuroimaging Initiative (ADNI) is an open-access dataset and can be accessed through their official website: [ADNI](https://adni.loni.usc.edu/).

- Citation:
  - **S. G. Mueller et al.**, “The Alzheimer’s Disease Neuroimaging Initiative,” *Neuroimaging Clinics of North America*, vol. 15, no. 4, p. 869, 2005.

## What is Included in This Repository
- All QPP templates from the paper, located in the `qpp_templates` folder:
  - **Stable Groups**:
    - sNC (Stable Normal Control)
    - sMCI (Stable Mild Cognitive Impairment)
    - sDAT (Stable Dementia of Alzheimer’s Type)
  - **Transitioning Groups**:
    - uNC-PRE (Unstable Normal Control - Pre-transition)
    - uNC-POST (Unstable Normal Control - Post-transition)
    - pMCI-PRE (Progressive Mild Cognitive Impairment - Pre-transition)
    - pMCI-POST (Progressive Mild Cognitive Impairment - Post-transition)
- Scripts for visualization and QPP template projection
- Statistical analysis outputs for Kruskal-Wallis and Dunn's post-hoc tests

## What You Need to Run the Code
1. **Preprocess the fMRI Scans**
   - The scans must first be preprocessed to remove artifacts and motion effects.
2. **Map Preprocessed Scans to an Atlas**
   - The NeuroMark2 atlas is used for this study ([Access NeuroMark2](https://trendscenter.org/data/)).
   - Citation:
     - **K. M. Jensen et al.**, “Addressing inconsistency in functional neuroimaging: A replicable data-driven multi-scale functional atlas for canonical brain networks,” *bioRxiv*, 2024.
     - **A. Iraji et al.**, “Identifying canonical and replicable multi-scale intrinsic connectivity networks in 100k+ resting-state fMRI datasets,” *Human Brain Mapping*, vol. 44, no. 17, pp. 5729–5748, 2023.
3. **Run QPPLab**
   - The preprocessed fMRI data must be processed through QPPLab, which extracts QPP templates and provides initial visualizations ([QPPLab Repository](https://github.com/GT-EmoryMINDlab/QPPLab)).
   - Citation:
     - **N. Xu et al.**, “QPPLab: A generally applicable software package for detecting, analyzing, and visualizing large-scale quasi-periodic spatiotemporal patterns (QPPs) of brain activity,” *bioRxiv*, 2023.

## How to Run the Code
1. **Run the Visualization Script**
   - This script will generate figures showing the different domains and subdomains of the QPPs included in the `qpp_templates` folder.
2. **Compute Correlations**
   - Users can compute the correlation of QPPs across groups to analyze functional connectivity differences.
3. **Project QPPs Onto Time-Series Data** (if available)
   - If time-series data is available, users can project the QPP templates onto their data to replicate the results.
4. **Perform Statistical Analysis**
   - The Kruskal-Wallis test is used to determine significant differences across groups.
   - Dunn’s test is used as a post-hoc analysis to identify significant group-level differences.

## Contact
For questions, please contact **Theodore J. LaGrow** at tlagrow@gatech.edu.

## License
This repository is shared under an open-access license for non-commercial use.
