# ShiftShareSEMatlab

Confidence intervals in shift-share designs (also called [Bartik
(1991)](http://research.upjohn.org/up_press/77/) designs) using procedures from
[Adão, Kolesár, and Morales (2018)](https://arxiv.org/abs/1806.07928). See the
[ShiftShareSE](https://github.com/kolesarm/ShiftShareSE) package for R version
of this code, and the
[ShiftShareSEStata](https://github.com/zhangxiang0822/ShiftShareSEStata) package
for a Stata version.

## Code description

### `ivreg_ss.m`

Implements the shift-share IV regression and reports AKM and AKM0 Confidence
Intervals.

```
[hat_beta, SE, pvalue, CIl, CIu, CItype] =
    ivreg_ss(Yn, Xn, Zn, controls, ln, weight,
                       cluster_vec, alpha, AKMtype, beta0)
```

#### Description of arguments

  `Yn`
  : dependent variable

  `Xn`
  : endogenous regressor

  `Zn`
  : shift-share IV

  `controls`
  : control matrix---vector of ones if empty

  `ln`
  : matrix of shares used in shift-share regressor. Rows are regions and columns
    are vectors. Notice that, for each row, the sum of columns must be less than
    or equal to 1 This matrix must have linearly independent columns.

  `weight`
  : observation weights

  `cluster_vec`
  : vector of clusters---no clustering if empty

  `alpha`
  : significance level for confidence interval

  `AKMtype`
  : 1 for AKM and 0 for AKM0

  `beta0`
  : null hypothesis

#### Description of output

  `hat_beta`
  : estimated coefficient on endogenous regressor

  `SE`
  : length of CI

  `pvalue`
  : p-value of the null hypothesis H0: beta = beta0

  `CIl`
  : lower bound of CI

  `CIu`
  : upper bound of CI

  `CI type`
  : 0 - AKM, 1 - standard AKM0, 2 - nonstandard AKM0 of the form
    `[-Inf,CIl]` U  `[CIu,Inf]`, 3 - nonstandard AKM0 of the form `[-Inf,Inf]`

### `reg_ss.m`

Implements the shift-share OLS regression and reports AKM and AKM0 Confidence
Intervals

```
[hat_beta, SE, pvalue, CIl, CIu, CItype] =
    reg_ss(Yn, Xn, controls, ln, weight,
                        sec_cluster_vec, alpha, AKMtype, beta0)
```

#### Description of arguments

All variables identical to `ivreg_ss.m`, except that `Xn` is the
shift-share regressor

#### Description of Output

Identical to `ivreg_ss.m`

### ADHapplication.m

Example to implement the code to generate the results in Column (2) of Table 5
in the paper. Uses the data `data_input_empADH.mat` and `endog_var_cz_ADH.mat`.
