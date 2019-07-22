Welcome to our summer research project website! This project focuses on studying ways to cluster longitudinal categorical data.

What is **longitudinal categorical data**?

- It consists of repeated measurements over time of an categorical variable, observed for many units or individuals. We'll call a set of measurements on one unit a *sequence*.  

Think about your sleep: every night’s sleep varies in length, but each moment can be classified into one of four categories: awake, REM sleep, light sleep and deep sleep. In this case, each unit is one night's sleep, which could be observed on the same individuals or on multiple individuals. 

It’s natural to think that the similarity between any two nights’ sleep or sequences can be quantified in many different ways. Our goal is to try and understand questions such as what factors determine the similarity between two sleep sequences and what sequences are similar enough to be considered in the same group or, more formally, “cluster”.

Alternatively, cluster analysis could be used to measure cancer patients’ end-of-life care. Starting the day of their diagnosis, cancer patients frequently transition between home, hospitals, skilled nursing facilities (SNF) and hospice. Given how patients transition between *spells*, we seek to identify patterns of healthcare utilization between patients that could help healthcare professionals make more informed decisions.

- The concept of *spells* is frequently used in sequence analysis. A spell is a longitudinal subsequence in which the subject stays in the same state. Take a sequence $A−A−B−B−B−A$, for example, each of the following are distinct spells: $A−A$, $B−B−B$, $A$.

Many methods are motiviated by genetics research. In this case, each position in the DNA sequence contains one of the four nitrogenous bases: *A*, *T*, *G*, *C*. Geneticists try to find the similarity between DNA sequences, but their cluster analysis methods can be extended to the realm of life trajectory analysis.

