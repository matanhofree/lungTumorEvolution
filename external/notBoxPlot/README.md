# notBoxPlot
<p align="center">
<img src="http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/26508/versions/17/screenshot.jpg"  height=350px/>
</p>


notBoxPlot is a MATLAB data visualisation function. 

Whilst box plots have their place, it's sometimes nicer to see all the data, rather than hiding them with summary statistics such as the inter-quartile range. This function (with a tongue in cheek name) addresses this problem. The use of the mean instead of the median and the SEM and SD instead of quartiles and whiskers are deliberate.
Jittered raw data are plotted for each group. Also shown are the mean, and 95% confidence intervals for the mean. This plotting style is designed to be used alongside parametric tests such as ANOVA and the t-test. Comparing the jittered data to the error bars provides a visual indication of whether the normality assumptions of the statistical tests are being violated. Furthermore, it allows one to eyeball the data to look for significant differences between means (non-overlapping confidence intervals indicate a significant difference at the chosen p-value, which here is 5%). Also see: http://jcb.rupress.org/cgi/content/abstract/177/1/7 Finally, 1 SD is also shown. Note that if data are not normally distributed then these statistics will be less meaningful.

The function has several examples and there are various visualization possibilities in addition to those shown in the above screenshot. For instance, the coloured areas can be replaced by lines.

Although it's worked well for situations I've needed it, I will be happy to modify the function if users come up against problems.

## Included functions
- notBoxPlot.m - generates plots as shown in screenshot
- SEM_calc.m - calculate standard error of the mean. Provided as a separate function file so that it can be used for other purposes.
- tInterval_Calc.m - calculate a t-interval. Right now notBoxPlot doesn't make use of this (unless the user edits the code, of course), but it still might be useful. For small sample sizes, the t-interval is larger than the SEM.

## Note
The MATLAB statistics toolbox is not required if you install the nantoolbox from here: http://pub.ist.ac.at/~schloegl/matlab/NaN/ Otherwise you will need the statistics toolbox for nan-handling.
