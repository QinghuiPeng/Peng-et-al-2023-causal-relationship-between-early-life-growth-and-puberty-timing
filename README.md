# Fetal genetically determined birth weight on puberty timing: a two-sample MR analysis
We used Mendelian randomization (MR) analyses to explore the causal relationships between fetal genetically determined birth weight (BW), childhood body mass index (BMI), and puberty timing (PT) and then explored the potential mediating role of childhood BMI.

## Environment details
R package to generate results as per Peng et al. 2023, Fetal genetically determined birth weight plays a causal role in earlier puberty timing: evidence from human genetic studies.

## Data sources
Data for BW and childhood BMI can be found at the following website: EGG Consortium http://egg-consortium.org/. Summary data for PT can be obtained from ReproGen Consortium https://reprogen.org/data_download.html.

## Analysis steps
1.Data pre-processing for BW, childhood BMI, and PT.

2.BW on PT & Childhood BMI on PT (Two-sample MR): Exposure: BW or childhood BMI; Outcome: PT.

3.Indirect effect (Two-step MR): Exposure: BW; Outcome: PT; Mediator: childhood BMI.

4.Direct effect (MVMR): Exposure: BW and childhood BMI; Outcome: PT.

5.Cluster (MR-Cluster): Exposure: BW; Outcome: PT.
