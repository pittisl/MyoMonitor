# MyoMonitor related source code

This repository holds related source code for the paper
"[MyoMonitor: Evaluating Muscle Fatigue with Commodity Smartphones](https://pittisl.github.io/publication/2021-myomonitor/)"
as published on [*Smart Health*](https://www.sciencedirect.com/journal/smart-health)
journal in 2021. You may find the paper [here](https://doi.org/10.1016/j.smhl.2020.100175).

## File Description

### Files under `MyoTester` directory

This is the project of the MyoTester Android App,
used to collect acoustic data using Android smartphone
with muscle workout.

* `App Usage.pdf`: The usage of Android App.

### Files under `muscle_v31` directory

These files process the collected data for fatigue
detection and classification.

* `Main.m`: Entry script. Classification of fatigue vs. non-fatigue on collected data.
* `Main2.m`: Classification of fatigue vs. non-fatigue on collected data, using classic approaches.
* `fatigueClassify.m`: the main fatigue classification algorithm.
* `featureExtraction.m`: Extract features from raw channel estimation.
* `preprocessing.m`: Preprocessing of received data by computing its channel estimation.
* `ReadAudioFile.m`: Script to read audio file to MATLAB data format.
* `sumDistanceFunc.m`: Function used by `preprocessing.m`, compute the sum distance of a complex data array to the centroid.
* `audioFeatureExtraction/` dir: third-party MATLAB library to provide some audio feature extraction functionalities.
* `accProcess.m`: Process accelerometer data collected during experiment.
* `accData/` dir: Sample input data to be processed by `accProcess.m`.
* `.mat` files: Workspace variables to be loaded.

### Files under `muscle_v1` directory

* `CircleFit.m`: Function used by `ModifiedMTI.m`.
* `curveSearch.m`: Function used by `singnalProcessMulti.m`.
* `ModifiedMTI.m`: Function used by `singalProcess.m`, `signalProcessBatch.m`, etc.
* `singalProcess.m` and related file: Process raw audio and plot channel estimation results.
* `Data/` dir: sample raw audio input for the channel estimation above.
* `sumDistanceFunc.m`: (duplicate file, see above)
* `ReadAudioFile.m`: (duplicate file, see above)
* `ORM.m`, `ORM_windowd.m`: Remove outliers in the fatigue data.
* `EMG/` dir: Sample of collected raw EMG signals.
* `transPointsIdentifyEMG.m`: Identify the transition point for EMG signal.
* `plotEMG.m`, `plotBothMulti.m`, and related file: plot collected EMG signals.
* `signalGeneration.m`: Generate audio signal for the experiment.
* `.mat` files: Workspace variables to be loaded by `signalProcess.m`.

### Other files

* `muscleInfoReadback_script.m`: Process PCM audio data collected from the phone.
* `kmeans_weight.m`: calculate k-means clustering for single muscle workout.
* `kmeans_segformat.m`: data segmentation for k means clustering.
* `channelEstiSigCreate.m`: create sequence to transmit 26-bit GSM training sequence.
* `dataSeg_timeseries.m`: Data segmentation for time series.
* `dataSeg.m`: Data segmentation for spectrum.

* * *

## License

Unless otherwise noted, all files under this repository are released with the following copyright information:

```
Copyright (c) 2020-2024 Intelligent Systems Lab, University of Pittsburgh. All Rights Reserved.
```

For files under `muscle_v31/audioFeatureExtraction/`, please refer to the individual license file within the directory.
