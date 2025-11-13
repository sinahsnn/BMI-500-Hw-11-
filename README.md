# title: BMI-500-Hw-11-
Model-based Machine Learning: Biopotential Modeling and Synthetic ECG Generation
---

# Author
**Name:** Mohammadsina Hassannia  
**Contact:** sina.hassannia@dbmi.emory.edu  
---

# Question Answered
**Question set 4** were selected. (Biopotential Modeling and Synthetic ECG Generation)

For full details, figures, and complete analysis, please refer to the accompanying PDF report:  
[Hassannia_BMI500HW11.pdf](Hassannia_BMI500HW11.pdf)
As an extra bonus, all the codes are modular and function based. 

---

# Key Insights and performance results

## Action Potential Dynamics & Oscillator Behavior (Part A)
- The second-order linear oscillator exhibits three regimes depending on damping parameter α:  
  **α < 0 → unstable diverging oscillations**,  
  **α = 0 →  Pure oscillatory behavior**,  
  **α > 0 → stable damped oscillations**.
  
- Phase-plane trajectories (outward spirals, closed ellipses, inward spirals) validate the theoretical stability predictions.
- **Extra:** Energy analysis shows exponential growth, conservation, or decay depending on α.

![Result A1](https://raw.githubusercontent.com/sinahsnn/BMI-500-Hw-11-/main/final_results/merged_a1.png)





## Self-Regulating Nonlinear Dynamics in the Van der Pol Oscillator (Part A)
- Nonlinear damping yields **stable limit-cycle oscillations**, controlling amplitude regardless of initial conditions.
- Increasing α moves behavior from sinusoidal → nonlinear → relaxation oscillations resembling action potentials.
- ω₀ scales oscillation frequency but not amplitude.
- Phase-plane changes demonstrate the transition from linear to highly nonlinear dynamics.
- **Extra:**: while the question asks for varying alpha, the variation of other parameter explored.

![Result A2](https://raw.githubusercontent.com/sinahsnn/BMI-500-Hw-11-/main/final_results/merged%20a2.png)


## FitzHugh–Nagumo vs Van der Pol (Part B)
- The FitzHugh–Nagumo model captures **fast–slow membrane dynamics**, excitability thresholds, and refractory periods.
- Increasing ε reduces time-scale separation, smoothing oscillations; varying I shifts behavior between quiescence and repetitive firing.
- The Van der Pol oscillator is simpler and produces similar relaxation oscillations, but FHN provides clearer physiological interpretation.
- - **Extra:**: while the question didnt ask for parameter variation variation of different parameters have been explored. for example it was explored that we should have at least I > I _threshold for firing.
![Result A3](https://raw.githubusercontent.com/sinahsnn/BMI-500-Hw-11-/main/final_results/merged%20a3.png)

## Synthetic ECG Modeling with McSharry–Clifford Model (Part B)
- The model produces physiologically realistic **P–QRS–T morphology** via limit-cycle oscillation and Gaussian shaping.
- Phase portraits (x–z, y–z, z–ż) show stable periodic loops representing cardiac cycles.
- **Extra:**: trajectories based on the paper has been explored.  

- Parameters **aᵢ** (amplitude) and **bᵢ** (width) modify morphology, enabling simulation of normal or pathological ECGs.
- - **Extra:**: the variable for all the components have been applied so that we proved the control on all the components while the question just ask for the effect of a and b.
- Heart-rate variability via ω(t) modulation creates realistic RR-interval fluctuations.
  ![result B](https://raw.githubusercontent.com/sinahsnn/BMI-500-Hw-11-/main/final_results/merged%20b.png)

## Stochastic Multichannel ECG Modeling & Clinical Realism (Part C)
- Stochastic VCG modeling introduces beat-to-beat variability and adds realistic noise sources:  
  **baseline wander**, **EMG muscle artifact**, and **Gaussian sensor noise**.  
- Noise processes mimic physiological and instrumentation disturbances seen in clinical ECGs.
- Resulting signals preserve morphology while showing natural variability, enhancing:
  - algorithm robustness testing  
  - signal-processing development  
  - machine-learning model training
- Stochastic modeling bridges idealized synthetic data with real-world clinical signals.
![Result C](https://raw.githubusercontent.com/sinahsnn/BMI-500-Hw-11-/main/final_results/sample%20c.png)

---

# Relevance to Model-Based Machine Learning

Model-based machine learning leverages mechanistic, physiologically driven models to get a better performance on the data.  The differential-equation systems explored—linear oscillators, Van der Pol, FitzHugh–Nagumo, and synthetic ECG generators—encode realistic temporal patterns and biophysical constraints.  These structured priors can reduce overfitting and help ML models generalize to real patient data.  Stochastic multichannel ECG modeling further enhances realism, allowing ML models to learn robustness to noise, artifacts, and variability, which are common challenges in clinical recordings. ( we can for example use ECG generation module to generate synthetic data if we do not have enough data. or we can make noisy data to make to model more robust) Overall, these models form a powerful foundation for training, benchmarking, and augmenting machine-learning systems in biomedical signal analysis.

---

# Suggestions for Future Modeling Improvements

- **Integrate the Hodgkin–Huxley model** to introduce detailed ionic-current dynamics and more biophysically grounded action potential modeling.
- **Develop dynamically coupled models based on cardiac dipole interactions** to simulate spatial conduction, propagation pathways, and coordinated atrial/ventricular activation.
- **Explore multi-scale modeling** that links HH ion-channel dynamics with dipole-level conduction for whole-heart electrophysiology.
- **Incorporate probabilistic or neural-ODE extensions** to blend mechanistic understanding with data-driven flexibility.
- **Expand stochastic modeling** to include electrode motion artifacts, nonstationary respiration, and adaptive noise processes.

** Note: Chatgpt has been used for some parts of this homework.(motly for checking my own codes and understanding the codes in the tollbox to be more clear, Part A.i, A.iii I used to checking and modeling  the problem. for part A.ii I used this for chcking. For part B, I did not use ChatGPT but I checked the codes with that. for part C I used OSET codes but I used GPT to summarize them and also understand how to use them. Also, as I want to implement similar code and just replace the functions I used chatGPT in this part.  
