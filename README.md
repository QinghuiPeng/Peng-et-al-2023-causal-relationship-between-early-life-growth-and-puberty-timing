# Peng-et-al-2023-causal-relationship-between-early-life-growth-and-puberty-timing
We used MR approaches to examine the causal relationships between birth weight (BW), childhood body mass index (BMI), and puberty timing (PT) and the potential mediating role of childhood BMI in the causal pathway.

## Environment details
R package to generate results as per Peng et al,. 2023, Deciphering the causal relationship between early-life growth and puberty timing: a Mendelian randomization study

## Data sources
Data for BW and childhood BMI can be found at the following website: EGG Consortium http://egg-consortium.org/. Summary data on PT can be obtained from ReproGen Consortium https://reprogen.org/data_download.html.

## Analysis steps
1.Data pre-processing for BW, childhood BMI, and PT.

2.BW on PT & Childhood BMI on PT (Two-sampe MR): exposure: BW or childhood BMI; outcome: PT.

3.Indirect effect (Two-step MR): exposure: BW; outcome: PT; mediator: childhood BMI.

4.Direct effect (MVMR): exposure: BW and childhood BMI; outcome: PT.

5.Cluster (MR-Cluster): exposure: BW; outcome: PT.
