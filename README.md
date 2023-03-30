
# GM-PHD WHISTLE DETECTOR

A Matlab set of functions and scripts for multi-target tracking of narrowband frequency modulated sounds- whistles. 

https://github.com/PinaGruden/GMPHD_whistle_contour_tracking

Copyright (c) 2016, Pina Gruden

The package provides a tool for multi-target tracking and extraction of dolphin whistle contours. The tracker is based on the Gaussian Mixture Probability Hypothesis Density Filter (GM-PHD) [Vo and Ma, 2006; Gruden and White, 2016].


## Author:
Pina Gruden - pina.gruden@gmail.com

## Technical references:

Vo, B.-N. and Ma, W.-K. (2006). The Gaussian mixture probability hypothesis density filter.
IEEE Transactions on Signal Processing, 54(11):4091-4104.

Gruden, P. and White., P. R. (2016). Automated tracking of dolphin whistles using Gaussian mixture probability hypothesis density filters. Journal of the Acoustical Society of America, 140(3):1981-1991.

## Modifications
---
I've introduced some modifications at the code. Actualy, now the code can be executed in cases when the file is much larger. You have to put a break point for visualize the results. Currently the script calculate the whistle candidates using the same contour whistle detection algorithm. The candidates are generated with two techniques:
* Picnogram
* Spectrogram
Both are introduced to the multi-tracking contour whistle algorithm and we save each one in an array of structures

## Important output variables and plots
---
We use a window of 1 sec for all the region analysis. On each one we can represent different plots:

- *Figure(5)*: Spectrogram of the time window.
- *Figure(4)*: Candidates and detections of spectrogram. 
- *Figure(3)*: Candidates and detections of picknogram. 

The usefull variables of the output are the next:
- **E_pic**: struct array with all the picknogram detections. They have two fiels:
    - time: time (in seconds) of the whistle
    - freq: frequencial points of the whistle (in Hz)
- **E_sp**: struct array with all the spectrogram detections. The same fields that the previous variable.


