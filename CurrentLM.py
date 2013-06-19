#
# current linear model is from midpointThresholdAnalysis.R
# lmTopTwoDecilesPlusC
#
#residuals shows some unbalanced tail effects, but mainly outside 1Q and 3Q
#Residual standard error: 6.921 on 439 degrees of freedom
#Multiple R-squared: 0.6907,	Adjusted R-squared: 0.6879 
#F-statistic: 245.1 on 4 and 439 DF,  p-value: < 2.2e-16
# test error metrics:
# MAE 6.1313; MEFR 2.2662
#
ileNames = ["8D", "9D"]
iles     = [.8, .9]
def applyCurrentLM(example, color):
    prediction = 16.8198 # intercept
    # B is intercept
    if color == "G":
        prediction += -12.0459
    elif color == "R":
        prediction += -15.2042
    prediction += example[color+"8D"] * -0.4423
    prediction += example[color+"9D"] *  0.8522
    return prediction
