# Classifying Plastic using mmWave Radar

This project aims to classify post-consumer plastic using commercial off-the-shelf (COTS) mmWave radar.

## How to use

This repository is configured for the AWR6843ISK.

First, set up the radar in MATLAB. MathWorks provides an interactive setup guide for supported boards here: https://www.mathworks.com/help/radar/ref/mmwaveradarsetup.html

Next, capture the beat/IF/IQ signal using `mmwave_capture.m`.

Then, visualize the captured data using `mmwave_viz_phase_ampl.m`. This file also includes a training script for a simple three-layer neural network to classify plastic types.

The dataset collected for my experiments is also included in this repository.
