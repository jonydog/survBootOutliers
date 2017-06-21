# survBootOutliers
This package provides three new outlier detection methods for survival datasets: BHT, OSD and DBHT.
The three methods propose to perform outlier detection in a
multivariate setting, using the Cox regression as the model and the concordance c-index as a measure of goodness
of fit. The first method BHT (Bootstrap Hypothesis Testing) is a sequential greedy procedure that presents a delete-1 statistic in order to make a bootstrap hypothesis test, testing for the increase in the concordance c-index. The second method OSD (One Step Deletion) is based on a sequential procedure that maximizes the c-index of the model using a greedy one-step-ahead search. The third method DBHT (Dual Bootstrap Hypothesis Testing) uses two antagonistic bootstrap resamples to perform the hypothesis test. These three
methods can be used to perform robust estimation for the Cox regression or other survival models, removing from the regression or other model, the most outlying observations.

# Related publications

Pinto J., Carvalho A. and Vinga S. (2015). Outlier Detection in Survival Analysis based on the Concordance C-index.In Proceedings of the International Conference on Bioinformatics Models, Methods and Algorithms - Volume 1: BIOINFORMATICS, (BIOSTEC 2015) ISBN 978-989-758-070-3, pages 75-82. DOI: 10.5220/0005225300750082

Pinto J.D., Carvalho A.M., Vinga S. (2015) Outlier Detection in Cox Proportional Hazards Models Based on the Concordance c-Index. In: Pardalos P., Pavone M., Farinella G., Cutello V. (eds) Machine Learning, Optimization, and Big Data. Lecture Notes in Computer Science, vol 9432. Springer, Chams
