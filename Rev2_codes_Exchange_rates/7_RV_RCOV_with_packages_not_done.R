rm(list=ls(all=TRUE))

# install.packages("realized", repos="http://R-Forge.R-project.org")

library(realized)

## Example

data(sbux.xts)
head(sbux.xts)
tail(sbux.xts)

data(lltc.xts)
head(lltc.xts)
tail(lltc.xts)

# is this log returns? or crude data? I think log returns, because many 0's

# https://rdrr.io/rforge/realized/man/rc.kernel.html
rc.kernel(x = sbux.xts,y = lltc.xts, kernel.param=1, kernel.type="bartlett", align.by ="seconds", align.period=60)

# https://rdrr.io/rforge/realized/man/rc.naive.html

rc.naive(x = sbux.xts, y=lltc.xts, period = 60, align.by ="seconds", align.period=1)

# https://rdrr.io/rforge/realized/man/rc.timescale.html

rc.timescale(x = sbux.xts, y=lltc.xts, period = 60,align.by ="seconds", align.period=1, adj.type="aa")






