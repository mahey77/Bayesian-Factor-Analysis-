library(ggplot2)

okun <- Okun_s_Law

okun <- na.omit(okun)



gokun.lm <- lm(gdp ~ ur, data = okun)
summary(okun.lm)
okun$yhat <- okun.lm$fitted.values
okun.predict <- predict.lm(okun.lm, interval = 'confidence', level = .95, se.fit = TRUE)
okun$upr <- okun.predict$fit[, 3]
okun$lwr <- okun.predict$fit[, 2]

g <- ggplot(okun, aes(y=gdp, x=ur)) + geom_point(col='blue') + ylim(-15,20)
g + xlab('Change in Unemployment Rate (percentage points)') + ylab('% Change in real GDP') + geom_line(aes(y=yhat), col='red') + geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.1)
