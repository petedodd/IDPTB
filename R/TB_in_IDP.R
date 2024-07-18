## === libraries
library(here)        #file references
library(refugees)    #IDP flow data
library(wpp2022)     #install_github("PPgp/wpp2022") : for crude death dates, CDR
library(countrycode) #converting country names
library(data.table)  #data handling
library(ggplot2)     #visualization
library(scales)      #visualization

## === utility functions
ssum <- function(x, ...) sqrt(sum(x^2, ...)) # for aggregating SDs (Pythagoras)
odds <- function(x) x / (1 - x)
iodds <- function(x) x / (1 + x)

## table formatting:
fmtB <- function(x, dg = 3, ...) formatC(signif(x, digits = dg), digits = dg, format = "fg", ...)
fmtS <- function(x) formatC(round(x, 1), digits = 2, format = "fg", flag = "#")
xyz <- function(x, y) paste0(fmtB(x), " (", fmtB(x - 1.96 * y), " to ", fmtB(x + 1.96 * y), ")")
xyw <- function(x, y) paste0(fmtS(x), " (", fmtS(x - 1.96 * y), " to ", fmtS(x + 1.96 * y), ")")



## === read in data
## crude death rates from WPP
data(misc1dt) # includes CDR
DD <- misc1dt[year == 2021, .(country_code, name, crdr = cdr / 1e3)]
DD[, iso3 := countrycode(name, origin = "country.name", destination = "iso3c")]
DD <- DD[!is.na(iso3)]


## UNHR data
D <- as.data.table(refugees::idmc)

## WHO TB estimates
fn <- here("data/TB_burden_countries_2023-11-08.csv")
T <- fread(fn)


## total TB incidence and mortality year yr w/ uncertainty
yr <- 2022 #set year of interest
totinc <- T[year == yr, .(
  M = sum(e_inc_num, na.rm = TRUE),
  S = ssum(e_inc_num_hi - e_inc_num_lo, na.rm = TRUE) / 3.92
)]
totmort <- T[year == yr, .(
  M = sum(e_mort_num, na.rm = TRUE),
  S = ssum(e_mort_num_hi - e_mort_num_lo, na.rm = TRUE) / 3.92
)]


## quick look at IDP
tot <- D[, .(tot = sum(total, na.rm = TRUE)), by = year]
ggplot(tot, aes(year, tot)) +
  geom_point() +
  geom_line() +
  theme_linedraw() +
  scale_y_continuous(label = comma, name = "Internally displaced total") +
  xlab("Year") +
  ggtitle("Last year = 2022")


ggsave(here("output/p.IDP.over.time.png"), h = 5, w = 7)


## ===  simple calcs/merging

## merge flows with TB data
BD <- merge(
  T[, .(year, iso3,
    ## baseline TB incidence and mortality per capita
    e_inc_100k, e_mort_100k,
    e_inc_100k.sd = (e_inc_100k_hi - e_inc_100k_lo) / 3.92,
    e_mort_100k.sd = (e_mort_100k_hi - e_mort_100k_lo) / 3.92,
    ## 2x incidnce
    e_inc_100kX2 = e_inc_100k * 2, e_mort_100kX2 = e_mort_100k * 2,
    ## 2x incidence AND 2x odds death
    e_inc_100kX22 = e_inc_100k * 2, e_mort_100kX22 = (e_inc_100k * 2) * iodds(2 * odds(cfr))
  )],
  D[year >= min(T$year), .(year, iso3 = coo_iso, coa_iso, total)], # flow data
  by = c("iso3", "year")
)
## merge in crude death rates
BD <- merge(BD, DD[, .(iso3, crdr)], by = "iso3", all.x = TRUE, all.y = FALSE)


## convert per capita rates to totals:
##    baseilne
BD[, tbi := e_inc_100k * total / 1e5]
BD[, tbm := e_mort_100k * total / 1e5]
BD[, tbi.sd := e_inc_100k.sd * total / 1e5]
BD[, tbm.sd := e_mort_100k.sd * total / 1e5]
##    2x incidence
BD[, tbi2 := e_inc_100kX2 * total / 1e5]
BD[, tbm2 := e_mort_100kX2 * total / 1e5]
##    2x incidence AND 2x OR death
BD[, tbi22 := e_inc_100kX22 * total / 1e5]
BD[, tbm22 := e_mort_100kX22 * total / 1e5]

## aggregate by year
BDSY <- BD[, .(
  tb = sum(tbi, na.rm = TRUE),
  tbd = sum(tbm, na.rm = TRUE),
  tb.sd = ssum(tbi.sd, na.rm = TRUE),
  tbd.sd = ssum(tbm.sd, na.rm = TRUE),
  tb2 = sum(tbi2, na.rm = TRUE),
  tbd2 = sum(tbm2, na.rm = TRUE),
  tb22 = sum(tbi22, na.rm = TRUE),
  tbd22 = sum(tbm22, na.rm = TRUE),
  totdeaths = sum(total * crdr * 2.5, na.rm = TRUE), # x factor 2.5 over background CDR
  total = sum(total)
), by = year]

## country
BC <- BD[year == yr]
BC <- merge(BC,
  T[
    year == yr,
    .(iso3, e_inc_num, e_inc_num.sd = (e_inc_num_hi - e_inc_num_lo) / 3.92)
  ],
  by = "iso3"
)
BC[, TBpc := 1e2 * tbi / e_inc_num]
BC[, TBpc.sd := TBpc * sqrt(tbi.sd^2 / tbi^2 + e_inc_num.sd^2 / e_inc_num^2)]
BCtop <- BC[rev(order(tbi))][1:10] #top 10 by incidence
BPCtop <- BC[rev(order(TBpc))][1:10] # top 10 by pc


## === trend estimation

## growth in movement
G <- BDSY[year >= 2018, .(total, year)]
G[, ratio := total / min(total)] # 10% / year for 2018-2022
mdl <- lm(total ~ year, data = G)
pp <- predict(mdl,
  newdata = data.frame(year = c(2025, 2030, 2035)),
  interval = "prediction"
) / G[year == 2022, total]
## ppS <- (pp[3]-pp[2])/3.92
pp <- as.data.table(pp)
## factor increases for 2025,2030,2035 vs 2022
pp[, S := (upr - lwr) / 3.92]


## TB trends
GTB <- T[, .(
  e_pop_num = sum(e_pop_num),
  e_inc_num = sum(e_inc_num),
  e_mort_num = sum(e_mort_num, na.rm = TRUE)
),
by = .(year)
]
GTB[, c("incpc", "mortpc") := .(e_inc_num / e_pop_num, e_mort_num / e_pop_num)]
GTBM <- melt(GTB, id = "year")


## NOTE COVID effects
ggplot(GTBM, aes(year, log(value), group = variable)) +
  geom_line() +
  facet_wrap(~variable, scales = "free") +
  stat_smooth(method = lm)



## regressions on each indicator
tlds <- list( # log transformed for safety
  pcinc = lm(log(incpc) ~ year, data = GTB),
  pcmort = lm(log(mortpc) ~ year, data = GTB),
  inc = lm(log(e_inc_num) ~ year, data = GTB),
  mort = lm(log(e_mort_num) ~ year, data = GTB)
)
tp <- list()
for (nm in names(tlds)) {
  tp0 <- predict(tlds[[nm]], newdata = data.frame(year = c(2025, 2030, 2035)), interval = "prediction")
  tp0 <- as.data.table(tp0)
  tp0[, c("fit", "lwr", "upr") := .(exp(fit), exp(lwr), exp(upr))]
  tp0[, qty := nm]
  tp0[, year := c(2025, 2030, 2035)]
  tp[[nm]] <- tp0
}
tp <- rbindlist(tp)

## ref values
refs <- GTB[year == 2020] # NOTE use 2020 as reference point due to COVID kinks
names(refs) <- c("year", "pop", "inc", "mort", "pcinc", "pcmort")
refs[, c("year", "pop") := NULL]
refs <- melt(refs)
## relatives
tp <- merge(tp, refs, by.x = "qty", by.y = "variable")
tp[, mid := fit / value]
## factor changes for TB indicators: 2025,2030,2035 vs 2020
tp[, S := (upr - lwr) / (3.92 * value)]


## === main table generation

tmp <- BDSY[year == yr] # data for 2022
## table data with projection
tmp2 <- rbindlist(
  list(
    ## baseline estimate 2022
    tmp[, .(
      assume = "baseline, 2022",
      TB = tb, TBM = tbd, TB.sd = tb.sd, TBM.sd = tbd.sd, totdeaths
    )],
    ## 2x incidence estimate 2022
    tmp[, .(
      assume = "2x incidence, 2022",
      TB = tb2, TBM = tbd2, TB.sd = 2 * tb.sd, TBM.sd = 2 * tbd.sd, totdeaths
    )],
    ## 2x incidence, 2x OR death estimate 2022
    tmp[, .(
      assume = "2x incidence, 2x OR death, 2022",
      TB = tb22, TBM = tbd22, TB.sd = 2 * tb.sd, TBM.sd = 4 * tbd.sd, totdeaths
    )],
    ## 2x incidence, 2x OR death estimate: projected 2025, 2030, 2035
    tmp[, .(
      assume = paste0("2x incidence, 2x OR death, ", c(2025, 2030, 2035)),
      ## (above) x (projected population) x (projected TB rates)
      TB = tb22 * pp[, fit] * tp[qty == "pcinc", mid],
      TBM = tbd22 * pp[, fit] * tp[qty == "pcmort", mid],
      ## uncertainties, using d(A*B)/(AB) = dlog(AB) = dA/A + dB/B etc
      TB.sd = tb22 * pp[, fit] * tp[qty == "pcmort", mid] *
        sqrt((2 * tb.sd)^2 / tb22^2 + pp[, S^2 / fit^2] + tp[qty == "pcinc", (S / mid)^2]),
      TBM.sd = tbd22 * pp[, fit] * tp[qty == "pcmort", mid] *
        sqrt((4 * tbd.sd / tbd22)^2 + pp[, S^2 / fit^2] + tp[qty == "pcmort", (S / mid)^2]),
      ## (IDP deaths) x (population increase)
      totdeaths = totdeaths * pp[, fit]
    )]
  )
)



## % of global TB totals
## incidence:
tmp2[, TBpc := 1e2 * TB / totinc$M]
tmp2[4:6, TBpc := 1e2 * TB / (totinc$M * tp[qty == "inc", mid])] #nums include trend; denominator trends
## deaths
tmp2[, TBMpc := 1e2 * TBM / totmort$M]
tmp2[4:6, TBMpc := 1e2 * TBM / (totmort$M * tp[qty == "mort", mid])] # denominator trend nums
## pc of all IDP deaths
tmp2[, pcall := 1e2 * TBM / totdeaths]
tmp2[4:6, pcall := 1e2 * TBM / (totdeaths * pp[, fit])]


## uncertainty
tmp2[, TBpc.sd := TBpc * sqrt((TB.sd / TB)^2 + (totinc$S / totinc$M)^2)]
tmp2[, TBMpc.sd := TBMpc * sqrt((TBM.sd / TBM)^2 + (totmort$S / totmort$M)^2)]
tmp2[, pcall.sd := 1e2 * TBM.sd / totdeaths]
## NOTE pp in num & denom: assume correlated so IDP growth not contributing unc'ty in fraction
tmp2[4:6, pcall.sd := pcall * sqrt((4 * tmp$tbd.sd / tmp$tbd22)^2 + tp[qty == "pcmort", (S / mid)^2])]
tmp2[, totdeaths := NULL] # remove again



## === write out results table:
table <- tmp2[, .(assume,
  `TB incidence` = xyz(TB, TB.sd),
  `TB deaths` = xyz(TBM, TBM.sd),
  `% of global TB incidence` = xyw(TBpc, TBpc.sd),
  `% of global TB deaths` = xyw(TBMpc, TBMpc.sd),
  `% of IDP deaths` = xyw(pcall, pcall.sd)
)]

table

fwrite(table, file = here("output/table.csv"))


## top 10 by incidence
top10 <- BCtop[, .(iso3, e_inc_100k,
                   IDP = total,
                   `TB incidence` = xyz(tbi, tbi.sd),
                   `TB deaths` = xyz(tbm, tbm.sd)
)]

fwrite(top10, file = here("output/top10.basecaseTB.csv"))

## top 10 by % of TB incidence
top10pc <- BPCtop[, .(iso3, e_inc_100k,
                      IDP = total,
                      `% of country TB incidence` = xyw(TBpc, TBpc.sd),
                      `IDP TB incidence` = xyz(tbi, tbi.sd),
                      `IDP TB deaths` = xyz(tbm, tbm.sd)
)]

fwrite(top10pc, file = here("output/top10pc.basecaseTB.csv"))
