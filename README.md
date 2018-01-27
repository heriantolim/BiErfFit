# BiErfFit
**BiErfFit** provides a tool to fit a hysteresis curve simultaneously with two linear combinations of error functions. For brevity, we call the curve which forms the left part of the hysteresis, *curve α*, and that which forms the right part, *curve β*. The model used to fit each curve is given by:

f(*x*) ~ `Baseline` + errorfunc1(*x*) + errorfunc2(*x*) + ... + errorfuncN(*x*)

where each error function is a derived form of the standard erf function as given by:

errorfunc(*x*) ~ (`Height`/2) (1 + erf( 2 sqrt(ln(2))(x - `Center`)/`Width`))

`Baseline`, `Center`, `Width`, and `Height` are the fit coefficients.

## Important Notes
- The algorithm in [fit.m](/@BiErfFit/fit.m) needs a major improvement, particularly for dealing with cases that involve curve fitting with more than one error functions.

## Licensing
This software is licensed under the GNU General Public License (version 3).

## Tested On
- MATLAB R2015b - R2017b

## Requirements
- MATLAB Curve Fitting Toolbox
- [PeakFit](https://github.com/heriantolim/PeakFit)

## Setting Up
1. Download or git-clone this repository and other repositories listed in the [Requirements](https://github.com/heriantolim/BiErfFit#requirements).
2. Add the repositories to the MATLAB's search path via `addpath(genpath( ... ))` OR this [version control system](https://github.com/heriantolim/MatlabVerCon).

## Usage
Construct the BiErfFit object in the following ways, and the fit results will be populated in the object's [public properties](https://github.com/heriantolim/BiErfFit#public-properties).

```MATLAB
obj=BiErfFit(Data, ... Name-Value ...)
```

OR

```MATLAB
obj=BiErfFit(XData, YData, ... Name-Value ...)
```

OR

```MATLAB
% Create an empty BiErfFit object.
obj=BiErfFit();

% Specify the data points and curve-fit settings via property assignments.
obj.XData= ... ;
obj.YData= ... ;
obj.Property1Name=Value1;
obj.Property2Name=Value2;
...

% Perform the curve fitting by reinstantiating the BiErfFit object.
obj=BiErfFit(obj);
```

`Data` must be specified as a cell containing two elements: the first for curve α and the second for curve β. Each element must be a two-column (or two-row) matrix where the first column (or first row) is the X data points and the second column (or second row) is the Y data points. In the alternative syntax, `XData` and `YData` are respectively the X and the Y data points, also specified as a two-element cell; but instead of matrix, the elements must be a vector.

The curve-fit settings can be specified after the mandatory arguments in Name-Value syntax, e.g. `BiErfFit(Data, 'Window', [10,90], 'NumErfs', [2,3])`. If no settings are specified, the algorithm will attempt to fit the hysteresis curve using the default settings. The default option fits curve α and curve β with one error function each, and may not produce a good fit. See [Best Practices](https://github.com/heriantolim/PeakFit#best-practices) for recommendations in making an optimal fit.

### Name-Value Pair Arguments
Any of the [public properties](https://github.com/heriantolim/BiErfFit#public-properties) can be specified as arguments in Name-Value syntax during the object construction. Specifying other things as Name-Value arguments will return an error.

## Examples

### Best Practices

## Public Properties
- `Data`, `XData`, `YData`: The data points of the curve to be fitted.
- `Window`: A vector of length two [*a*,*b*] that limits the fitting to only the data points whose X coordinates lies within [*a*,*b*].
- `NumErfs`: The number of error functions wished to be fitted for curve α and curve β. Specified as a positive integer vector of length two, [*n<sub>α</sub>*,*n<sub>β</sub>*]. When the corresponding start points, lower, or upper bounds are set with vectors of length greater than *n<sub>α</sub>* (or *n<sub>β</sub>*), then *n<sub>α</sub>* (or *n<sub>β</sub>*) will be incremented to adjust to the maximum length of these vectors. When the maximum length of these vectors is less, then these vectors will be expanded and filled with the default values. `NumErfs` defaults to [1,1].

### Fit Results
- (Read-only) `Center`, `Height`, `Width`: A cell containing two elements: a 3-by-*n<sub>α</sub>* matrix and a 3-by-*n<sub>β</sub>* matrix. Each matrix stores the fit results for the center, height, and width, respectively, for curve α or curve β. Each column in the matrices correspond to the fit results for a particular error function. The first row holds the values at convergence, and the second (third) row holds the 95% CI lower (upper) bounds.
- (Read-only) `Baseline`: A 3-by-1 vector that stores the fit results for the baseline. The first row holds the values at convergence, and the second (third) row holds the 95% CI lower (upper) bounds.
- (Read-only) `Asymptote`: The lower and upper asymptotes to which the combined error functions converge at *x* = ±∞. Given as a 3-by-2 matrix, where the first (second) column holds the values for the lower (upper) asymptote, the first row holds the values at convergence, and the second (third) row holds the 95% CI lower (upper) bounds.
- (Read-only) `DynamicRange`: The ratio of the increase or decrease in *y* as *x* goes from -∞ to +∞. Given as a 3-by-1 vector, where the first row holds the values at convergence, and the second (third) row holds the 95% CI lower (upper) bounds.
- (Read-only) `HysteresisWidth`: The hysteresis width. Given as a 3-by-1 vector, where the first row holds the values at convergence, and the second (third) row holds the 95% CI lower (upper) bounds.
- (Read-only) `TransitionWidth`: The transition width. Given as a cell containing two elements. Each element is a 3-by-1 vector holding the values for curve α or curve β, respectively. The first row in each vector holds the values at convergence, and the second (third) row holds the 95% CI lower (upper) bounds.
- (Read-only) `RelStDev`: The relative standard deviation of the fit results.
- (Read-only) `CoeffDeterm`: The coefficient of determination of the fit results.
- (Read-only) `AdjCoeffDeterm`: The degree-of-freedom adjusted coefficient of determination of the fit results.
- (Read-only) `NumFunEvals`: The number of function evaluations.
- (Read-only) `NumIters`: The number of iterations.
- (Read-only) `ExitFlag`: Describes the exit condition of the algorithm. Positive flags indicate convergence, within tolerances. Zero flags indicate that the maximum number of function evaluations or iterations was exceeded. Negative flags indicate that the algorithm did not converge to a solution.

### Fitting Start Points
- `CenterStart`, `WidthStart`, `HeightStart`: The start points for the center, width, and height, respectively. Specified as a cell containing two vectors: the first for curve α and the second for curve β. Can also be specified as a vector of length two, if `NumErfs` is intended to be [1,1].
- `BaselineStart`: The initial value for the baseline coefficient. Specified as a real scalar.

### Fitting Lower Bounds
- `CenterLow`, `WidthLow`, `HeightLow`: The lower bounds for the center, width, and height, respectively. Specified as a cell containing two vectors: the first for curve α and the second for curve β. Can also be specified as a vector of length two, if `NumErfs` is intended to be [1,1].
- `BaselineStart`: The lower bound for the baseline coefficient. Specified as a real scalar.

### Fitting Upper Bounds
- `CenterLow`, `WidthLow`, `HeightLow`: The upper bounds for the center, width, and height, respectively. Specified as a cell containing two vectors: the first for curve α and the second for curve β. Can also be specified as a vector of length two, if `NumErfs` is intended to be [1,1].
- `BaselineStart`: The upper bound for the baseline coefficient. Specified as a real scalar.

### Algorithm Parameters
- (Read-only) `Method` = 'NonLinearLeastSquares'; The method used for the fitting.
- `Robust`: The type of the least-squares method to be used in the fitting. Avaliable values are 'off', 'LAR' (least absolute residual method), and 'Bisquare' (bisquare weight method). Defaults to 'off'.
- `Algorithm`: The algorithm to be used in the fitting. Available values are 'Lavenberg-Marquardt', 'Gauss-Newton', or 'Trust-Region'. Defaults to 'Trust-Region'.
- `MovMeanWidth`: The window width of the moving average used for smoothing the curve in order to filter out the noise before finding the maximas. This parameter is used only when CenterStart is not given. The value of this can be set as a positive integer which specifies the width in terms of the number of data points, OR a real scalar between [0,1] which specifies the width as a fraction of the total number of data points. Defaults to 0.02.
- `DiffMaxChange`: The maximum change in coefficients for finite difference gradients. Defaults to 0.1.
- `DiffMinChange`: The minimum change in coefficients for finite difference gradients. Defaults to 1e-8.
- `MaxFunEvals`: The allowed maximum number of evaluations of the model. Defaults to 1e5.
- `MaxIters`: The maximum number of iterations allowed for the fit. Defaults to 1e3.
- `TolFun`: The termination tolerance on the model value. Defaults to 1e-6.
- `TolX`: The termination tolerance on the coefficient values. Defaults to 1e-6.

## Public Methods
- `disp`: Displays the options, results, error and performance of the fitting.
- `model`: Returns the reconstructed data points (model) using the fit results. See [@BiErfFit/model.m](/@BiErfFit/model.m) for more info.

### Static Methods
- `fnerf`: The error function. See [@BiErfFit/fnerf.m](/@BiErfFit/fnerf.m) for more info.
- `fngaussian`: The Gaussian function. See [@BiErfFit/fngaussian.m](/@BiErfFit/fngaussian.m) for more info.
