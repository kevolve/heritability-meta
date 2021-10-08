---
title: "Coral Heritability Meta-analysis"
subtitle: "Supplementary Code Documentation B: Model Selection and Diagnostics"
author: "Kevin R Bairos-Novak"
date: "Document last run on 2021-10-08"
output:
  html_document: 
    df_print: kable
    toc: yes
    highlight: tango
    theme: cerulean
    number_sections: no
    fig_width: 10
    fig_height: 8
    fig_caption: yes
    code_folding: hide
    self_contained: TRUE
    keep_md: yes
editor_options: 
  chunk_output_type: console
knit: (function(inputFile, encoding) { rmarkdown::render(inputFile, encoding = encoding, output_file = paste0(xfun::sans_ext(inputFile),"_",format(Sys.time(), '%Y-%m-%d'), ".html")) })
---



```{=html}
<style>
.container { width: 1000px; }
h2 { color: #f8f8f8; background-color: #437FAA; }
h3 { color: #f8f8f8; background-color: #437FAA; text-align: center; }
<!-- Get highlight from Terminal: pandoc --print-highlight-style tango -->
</style>
```


# Start


```{.r .fold-show}
# rm(list=ls())
setwd("~/Documents/PhD Thesis/Heritability meta-analysis")

library(readxl)
library(tidyverse); theme_set(theme_light())
library(metafor)
library(patchwork)
source('Functions/useful_functions.R')
source('Functions/mlm.variance.distribution.R')

# Load processed data:
load("Data/heritability_estimates_processed.RData")
# data.full <- data.full %>% filter(!(study == "Zhang et al.")) # duplicate estimate!
data <- data.full %>% 
	filter(trait != "gamete contribution") %>% 
	mutate(trait = fct_drop(trait))
dat_s <- data %>% 
	filter(!(trait %in% c("photochemistry", "immune response"))) %>% 
	mutate(trait = factor(trait))
dat_sh <- data %>% 
	filter(!(trait %in% c("photochemistry", "nutrient content", "immune response"))) %>% 
	mutate(trait = factor(trait))
dat_shg <- data %>% 
	filter(trait %in% c("symbiont community", "survival", "nutrient content", "growth", "bleaching")) %>% 
	filter(!(growth.form %in% c("columnar", "encrusting"))) %>%
	mutate(trait = factor(trait), growth.form = factor(growth.form))
dat_tm <- data %>% 
	filter(!is.na(temp.diff)) %>%
	filter(!(trait %in% c("symbiont community", "morphology", "gene expression"))) %>% 
	mutate(trait = factor(trait),
		   temp.manip = temp.manip.plus.others,
		   temp.s = sqrt(temp.diff))
dat_tm_nz <- dat_tm %>% filter(temp.diff >0) %>%
	mutate(trait = factor(trait), 
		   stage = factor(stage), 
		   h2 = factor(h2),
		   growth.form = factor(growth.form))


# Specify final main model
final.model.A <- rma(yi=val.log ~ trait*stage + h2, 
					 vi=sv.log, data=dat_s, method="REML", test="knha")
```

------------------------------------------------------------------------

# Summary of results

When defining the base model architecture that includes all data (save
for one estimate on the heritability of gamete compatibility), the
optimal and most parsimonious explanatory variable was the trait type.
In order of increasing relative heritability (not accounting for
differences in heritability type):

-   gene expression
-   photochemistry
-   growth
-   nutrient content
-   bleaching
-   morphology
-   symbiont community
-   immune response
-   survival

In sub-analyses using smaller data subsets, life stage, heritability type, and growth form
came out important, whereas temperature manipulation
(binary yes/no), temperature differential (°C above
normal/ambient/control temperature), and symbiont type (symbiotic or
asymbiotic) were not.

An overall model of trait type, and a sub-analysis model of trait type x life stage + heritability type was preferred through model selection. No other models were found. 


# Overview of sample sizes

Narrow-sense heritability is reported for 10 studies, broad-sense for 10
studies, including one study reporting broad-sense and narrow-sense
estimates (Carlon et al. 2011).


```r
table.plot("h2", "study", data.full) # only Carlon et al. study with multiple h2 values
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


```r
summarise.coverage("trait", "study", data.full)
```

Most studies (12/19) report on a single trait type and not multiple
trait types. Three studies report on two trait types (Meyer et al. 2009;
Manzello et al. 2019; Kenkel et al. 2015), two studies report on three
trait types (Csázár et al. 2010; Quigley et al. 2020), and Zhang et al.
(2019) and Wright et al. (2019) report on four and six trait types,
respectively.

Both immune response and gamete contribution trait types have only a
single study contributing to their estimates. In the case of the latter,
gamete contribution has only one estimate, limiting any inferences using
models.


```r
table.plot("trait", "species", data.full)
data.full %>% group_by(species) %>% tally %>% arrange(-n)
summarise.coverage("species", "study", data.full)$Y
```

Most heritability estimates are from Acropora millepora (48 estimates; 5
separate studies), A. spathulata (12; 1), Porites astreoides (9; 2),
Fabia fragum (6; 1), Orbicella faveolata (4; 3), and A. cervicornis (3;
2).

Most diversity of species/studies are covered by trait types such as
growth, survival, bleaching, and symbiont community.

# Model Selection Process

## Complete Dataset

Following the process of mixed model selection outlined by Zuur et al.
(Mixed Effects Modelling Book, 2009, pp. 121-122) and Diggle et al.
(2002), we:

1.  Begin with a 'beyond optimal' model with as many fixed effects as
    possible, given the available data
2.  With the beyond optimal fixed effects, test for the optimal random
    effects structure by comparing AIC among models fit using REML
3.  Using the optimal random effects structure, compare fixed effects by
    comparing AIC among models fit using ML
4.  Refit the final mixed effects model using REML
5.  Check final model diagnostics: funnel plots, fail safe number, Cook's distances, residual plots, and simulated data

### Step 1: Begin with a *beyond optimal* model

We have 4 main factors that cover the complete dataset:

-   trait type (9 levels, including survival, growth, bleaching, etc.)

-   heritability type (broad- or narrow-sense)

-   coral life stage (larval, juvenile, adult)

-   coral growth form (3 main types: corymbose, branching, massive; also
    encrusting and columnar)

Thus, our beyond optimal model is one that includes the main effects of
all 4 factors (additive model). As seen later on, we cannot examine
interactions, given factor combinations being missing For example,
heritability estimates for narrow-sense heritability exist only for 6/9
trait types, so no interaction between trait type x heritability type
can be properly examined without removing the 3 trait types with only
broad-sense heritability estimates.


```r
summarise.coverage("trait", "h2", data)$X # 9 broad, 6 narrow
```

### Step 2: Select random effects structure

Next, we fit different random effect structures, including nested
(3-level) models that include either study ID or species ID as an upper
level, with estimate ID nested within (sampling variance acts as the
3rd, most basal level in a meta-model). We also fit similar models with
one or both of study/species ID and estimate ID fixed at zero, as well
as a 'null' random effects model, where all random effects are fit at
zero (essentially a fixed effects model).


```r
m.study.nested <- rma.mv(yi=val.log ~ trait + h2 + growth.form + stage,
						 V=sv.log, random = ~1|study/est.id, 
						 data=data, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.spp.nested <- rma.mv(yi=val.log  ~ trait + h2 + growth.form + stage,
						 V=sv.log, random = ~1|species/est.id, 
						 data=data, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.only.est <- rma.mv(yi=val.log ~ trait + h2 + growth.form + stage,
						 V=sv.log, random = ~1|study/est.id, 
						 data=data, method="REML", tdist=T,
						 sigma2 = c(0, NA)) # same as random = ~1|est.id
m.only.study <- rma.mv(yi=val.log ~ trait + h2 + growth.form + stage,
					  V=sv.log, random = ~1|study/est.id, 
					  data=data, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|study
m.only.spp <- rma.mv(yi=val.log ~ trait + h2 + growth.form + stage,
					  V=sv.log, random = ~1|species/est.id, 
					  data=data, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|species
m.no.reff <- rma.mv(yi=val.log ~ trait + h2 + growth.form + stage,
					  V=sv.log, random = ~1|study/est.id, 
					  data=data, method="REML", tdist=T,
					  sigma2 = c(0, 0)) # no random effects

TableS1 <- myAIC(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff) %>%
	mutate(model = fct_recode(model,
		"Estimate ID nested in study ID" = "m.study.nested",
		"Estimate ID nested in species" = "m.spp.nested",
		"Study ID only" = "m.only.study",
		"Species only" = "m.only.spp",
		"Estimate ID only" = "m.only.est",
		"No random effect" = "m.no.reff"))
write.csv(TableS1, file="Suppl Tables//TableS1.csv", row.names=F)
TableS1
```

<div class="kable-table">

|model                          | df|     AICc|     ΔAICc|
|:------------------------------|--:|--------:|---------:|
|Estimate ID nested in study ID | 15| 112.2116|   0.00000|
|Estimate ID only               | 15| 118.1781|   5.96650|
|Estimate ID nested in species  | 15| 121.5713|   9.35972|
|Study ID only                  | 15| 190.3282|  78.11664|
|Species only                   | 15| 311.9348| 199.72324|
|No random effect               | 15| 353.5841| 241.37256|

</div>

Here, we see that the model with estimate ID nested within study ID is
preferred (`random=~1|study/est.id, sigma2 = c(NA,NA)`), since it has
the lowest relative AICc value.

Thus, we continue on with this random effects structure when deciding
our fixed effects structure.

### Step 3: Select fixed effects structure

Next, we select the optimal by comparing model AICc values of all
possible fixed-effect models for factors that exist for the complete
dataset:

-   trait type (*t*; 9 levels: survival, growth, gene expression, etc.)

-   heritability type (*h*; 2 levels: $h^2$ or $H^2$)

-   coral growth form (*g*; 4 levels: branching, massive/columnar,
    corymbose, or encrusting)

-   coral life stage (*s*; 3 levels: larval, juvenile, adult)

We also have a factor of symbiont type (symbiotic or asymbiotic),
however, there are only 4 estimates with asymbiotic individuals, and
these all involve coral larvae from broadcast-spawners. Thus, there are
too few estimates and too little coverage to compare at this stage.
Similarly, temperature was not manipulated for a few of the studies, and
thus is not assessed here either.

We fit all combinations of additive factors possible, excluding any
interactions at this point, as the data are not complete for a number of
factor combinations. We display only the top 10 models by AICc (given it
is different in all models to AIC, our N sample size is too small for
the k model parameters being estimated!).

Again, because the data are not complete for interactions at this stage,
we examine only additive models with the various factors.


```r
dat <- data %>% mutate(t = trait, h = h2, g = growth.form, s = stage)
mod_names <- ls(pattern = "m.fixed.")
do.call("rm", as.list(mod_names))

# One-ways with less and less other terms:
m.fixed.1way.t.h.g.s <- rma.mv(yi=val.log ~ t + h + g + s, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.t.h.g <- rma.mv(yi=val.log ~ t + h + g, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.t.h.s <- rma.mv(yi=val.log ~ t + h + s, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.t.g.s <- rma.mv(yi=val.log ~ t + g + s, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.h.g.s <- rma.mv(yi=val.log ~ h + g + s, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.t.h <- rma.mv(yi=val.log ~ t + h, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.t.g <- rma.mv(yi=val.log ~ t + g, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.t.s <- rma.mv(yi=val.log ~ t + s, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.h.g <- rma.mv(yi=val.log ~ h + g, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.h.s <- rma.mv(yi=val.log ~ h + s, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.g.s <- rma.mv(yi=val.log ~ g + s, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)

m.fixed.1way.t <- rma.mv(yi=val.log ~ t, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.h <- rma.mv(yi=val.log ~ h, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.g <- rma.mv(yi=val.log ~ g, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.s <- rma.mv(yi=val.log ~ s, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)
m.fixed.1way.null <- rma.mv(yi=val.log ~ 1, 
								 random = ~1|study/est.id, V=sv.log,
								 data=dat, method="ML", tdist=T)


mod_names <- ls(pattern = "m.fixed")
param_list=c("t", "h", "g", "s")
table1 <- do.call(myAIC, as.list(c(mod_names, getmodel=T, 
									 param.list = list(param_list))))
# table1$rel.importance

AICt <- table1$AICtable %>% filter(model == "m.fixed.1way.t") %>% pull(4) %>% round(2)
AICth <- table1$AICtable %>% filter(model == "m.fixed.1way.t.h") %>% pull(4) %>% round(2)
AICts <- table1$AICtable %>% filter(model == "m.fixed.1way.t.s") %>% pull(4) %>% round(2)
AICths <- table1$AICtable %>% filter(model == "m.fixed.1way.t.h.s") %>% pull(4) %>% round(2)

# head(table1$AICtable)

TableS2 <- head(table1$AICtable) %>%
	mutate(model = fct_recode(model,
		"trait" = "m.fixed.1way.t",
		"trait + heritability type" = "m.fixed.1way.t.h",
		"trait + life stage" = "m.fixed.1way.t.s",
		"trait + heritability type + life stage" = "m.fixed.1way.t.h.s",
		"trait + growth form + life stage" = "m.fixed.1way.t.g.s",
		"trait + growth form" = "m.fixed.1way.t.g"))
write.csv(TableS2, file="Suppl Tables/TableS2.csv", row.names=F)
TableS2
```

<div class="kable-table">

|model                                  | df|     AICc|     ΔAICc|
|:--------------------------------------|--:|--------:|---------:|
|trait                                  |  8| 100.7705| 0.0000000|
|trait + heritability type              |  9| 101.1370| 0.3664769|
|trait + life stage                     | 10| 101.8590| 1.0884955|
|trait + heritability type + life stage | 11| 103.3778| 2.6072277|
|trait + growth form                    | 12| 109.0199| 8.2493425|
|m.fixed.1way.t.h.g                     | 13| 109.3342| 8.5636276|

</div>

The model with the lowest AICc is one with only trait type as a fixed
effect, followed closely (within ΔAICc \< 2) by a model of trait + life stage (ΔAICc =
1.09), and trait +
heritability type (ΔAICc = 0.37), and finally, trait + heritability + stage (ΔAICc =
2.61). and finally a few models including life stage,
heritability, and trait type. However, clearly, there are data
limitations present that limit the number of model parameters possible.
Thus, we accept trait type by itself as the core model for the complete
dataset.


```r
rbind(
	data.frame(
		Model = rep(paste0("~ trait (ΔAICc=", AICt,")")),
		Parameter = names(coef(m.fixed.1way.t)), 
		Estimate = coef(m.fixed.1way.t),
		SE = m.fixed.1way.t$se),
	data.frame(
		Model = rep(paste0("~ trait + life stage (ΔAICc=", AICts,")")),
		Parameter = names(coef(m.fixed.1way.t.s)), 
		Estimate = coef(m.fixed.1way.t.s),
		SE = m.fixed.1way.t.s$se),
	data.frame(
		Model = rep(paste0("~ trait + h2 (ΔAICc=", AICth,")")),
		Parameter = names(coef(m.fixed.1way.t.h)), 
		Estimate = coef(m.fixed.1way.t.h),
		SE = m.fixed.1way.t.h$se),
	data.frame(
		Model = rep(paste0("~ trait + h2 + life stage (ΔAICc=", AICths,")")),
		Parameter = names(coef(m.fixed.1way.t.h.s)), 
		Estimate = coef(m.fixed.1way.t.h.s),
		SE = m.fixed.1way.t.h.s$se)) %>%
	mutate(Parameter = factor(Parameter),
		   Parameter = fct_relevel(Parameter, "intrcpt")) %>%
	ggplot(aes(x=Estimate, y=Parameter, color=Model)) +
	geom_vline(xintercept=0, linetype="dashed") +
	geom_pointrange(aes(xmin = Estimate - SE, xmax = Estimate + SE), position = position_dodge(width=0.3)) + 
	labs(x="Parameter estimate (log[X+0.2] scale) ± SE") +
	theme(legend.position="top", legend.direction="vertical")
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

The plot shows relatively similar coefficient estimates and uncertainty
in the top 4 models. However, we go with the best model of only trait
type as a predictor.

### Step 4: Fit final model with REML

Due to the small effect size of life stage in the trait + life stage model, and the nearly identical AICc value (difference of 0.15), we fit our final model using trait type only as:

$$ logit(h^2 + 0.2) \sim N(\mu, \sigma^2) \\
\mu = trait_k \,  + \, \epsilon_{i} $$

where $\sigma^2$ are the estimated amount of heterogeneity in the
studies/estimates, calculated using an inverse variance weighting method
via package `metafor` that accounts for uncertainty in each estimate,
and $\epsilon$ is the residual error for each estimate.


```r
# Using subset data
final.model <- rma.mv(yi=val.log ~ trait,
						 V=sv.log, random = ~1|study/est.id, 
						 data=data, method="REML", tdist=T)
final.model
```

```

Multivariate Meta-Analysis Model (k = 94; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0933  0.3055     18     no         study 
sigma^2.2  0.0566  0.2379     94     no  study/est.id 

Test for Residual Heterogeneity:
QE(df = 85) = 477.5080, p-val < .0001

Test of Moderators (coefficients 2:9):
F(df1 = 8, df2 = 85) = 6.5312, p-val < .0001

Model Results:

                         estimate      se     tval    pval    ci.lb    ci.ub 
intrcpt                   -1.1454  0.2175  -5.2667  <.0001  -1.5779  -0.7130 
traitphotochemistry        0.3631  0.1832   1.9818  0.0507  -0.0012   0.7273 
traitgrowth                0.3884  0.2140   1.8149  0.0731  -0.0371   0.8139 
traitnutrient content      0.5335  0.2375   2.2466  0.0273   0.0614   1.0056 
traitbleaching             0.5690  0.2438   2.3340  0.0220   0.0843   1.0537 
traitsymbiont community    0.6742  0.3090   2.1818  0.0319   0.0598   1.2886 
traitmorphology            0.7919  0.3083   2.5683  0.0120   0.1789   1.4050 
traitimmune response       0.9435  0.2532   3.7265  0.0003   0.4401   1.4468 
traitsurvival              1.1345  0.2199   5.1598  <.0001   0.6973   1.5716 
 
intrcpt                  *** 
traitphotochemistry        . 
traitgrowth                . 
traitnutrient content      * 
traitbleaching             * 
traitsymbiont community    * 
traitmorphology            * 
traitimmune response     *** 
traitsurvival            *** 

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
save(final.model, file="Models/final.model.Rdata")

x <- mlm.variance.distribution(final.model, suppress.plot=F)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

------------------------------------------------------------------------

### Step 5: Check final model diagnostics

Lastly, we check the final model diagnostics, such as examining the
funnel plot for asymmetry indicating publishing bias, in addition to simulating data
based on the final model to examine if it loosely resembles the data
observed.

#### Funnel plot


```r
x <- resid(final.model)
hist(x) # looks approx. normal!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
## Model checks:
par(mfrow=c(2,1))
funnel(final.model)
funnel(final.model, yaxis="vinv")
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

```r
par(mfrow=c(1,1))

ranktest(final.model) # p = 0.72
```

```

Rank Correlation Test for Funnel Plot Asymmetry

Kendall's tau = -0.0252, p = 0.7194
```

Visually, the funnel plot appears to be slightly asymmetric in favour of studies reporting lower heritabilities, but a rank correlation test for funnel plot asymmetry was not statistically significant (Kendall's $\tau = = -0.025$, $P = 0.72$), indicating no apparent publication bias.

#### Fail-safe number


```r
fsn(yi=val.log, vi=sv.log, data=data, type="Rosenberg") # N=1287 to get no significance (but assumes fixed effects)
```

```

Fail-safe N Calculation Using the Rosenberg Approach

Average Effect Size:        -0.0875
Observed Significance Level: <.0001
Target Significance Level:   0.05

Fail-safe N: 1285
```

```r
# Number of studies x5 + 10:

length(unique(data$study))*5 + 10
```

```
[1] 100
```

#### Cook's distance


```r
cdist <- cooks.distance(final.model, reestimate=FALSE)
cdist.i <- which(cdist >2)
plot(cdist)
points(cdist.i, cdist[cdist.i], col="red", pch=16)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
data[cdist.i,] %>% select(study, species, val, trait) %>%
	cbind(., cooks = round(cdist[cdist.i], 2))
```

<div class="kable-table">

|   |study              |species            |  val|trait           | cooks|
|:--|:------------------|:------------------|----:|:---------------|-----:|
|16 |Wright et al. 2019 |Acropora millepora | 0.92|immune response |   5.2|

</div>

```r
update(final.model, .~., data[-cdist.i,])
```

```

Multivariate Meta-Analysis Model (k = 94; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0933  0.3055     18     no         study 
sigma^2.2  0.0566  0.2379     94     no  study/est.id 

Test for Residual Heterogeneity:
QE(df = 85) = 477.5080, p-val < .0001

Test of Moderators (coefficients 2:9):
F(df1 = 8, df2 = 85) = 6.5312, p-val < .0001

Model Results:

                         estimate      se     tval    pval    ci.lb    ci.ub 
intrcpt                   -1.1454  0.2175  -5.2667  <.0001  -1.5779  -0.7130 
traitphotochemistry        0.3631  0.1832   1.9818  0.0507  -0.0012   0.7273 
traitgrowth                0.3884  0.2140   1.8149  0.0731  -0.0371   0.8139 
traitnutrient content      0.5335  0.2375   2.2466  0.0273   0.0614   1.0056 
traitbleaching             0.5690  0.2438   2.3340  0.0220   0.0843   1.0537 
traitsymbiont community    0.6742  0.3090   2.1818  0.0319   0.0598   1.2886 
traitmorphology            0.7919  0.3083   2.5683  0.0120   0.1789   1.4050 
traitimmune response       0.9435  0.2532   3.7265  0.0003   0.4401   1.4468 
traitsurvival              1.1345  0.2199   5.1598  <.0001   0.6973   1.5716 
 
intrcpt                  *** 
traitphotochemistry        . 
traitgrowth                . 
traitnutrient content      * 
traitbleaching             * 
traitsymbiont community    * 
traitmorphology            * 
traitimmune response     *** 
traitsurvival            *** 

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The main influential point is one for immune response, which drives the
average estimate for immune response much higher than other values for
immune response. Note that all the values for the immune response trait
come from the same study (Wright et al. 2019), and so the values are
already to be interpreted cautiously.

#### Residuals and simulated data


```r
plot(fitted(final.model), rstandard(final.model)$z)
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
plot(residuals(final.model, type="pearson"))
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-14-2.png)<!-- -->

```r
# Simulate new data given the parameter values of the data and same variance structure, plot on real data to see if it is similar
xsim <- simulate(final.model, 100) %>%
	mutate(true.val = data$val.log,
		   weighting = 1/data$sv.log) %>%
	pivot_longer(cols=-c(true.val, weighting)) %>%
	select(name, true.val, sim.val = value, weighting)
xsim %>%
	ggplot(aes(x=true.val, y=sim.val)) +
	geom_point(aes(col=weighting), alpha=0.05) +
	geom_abline(slope=1, intercept=0, linetype="dashed") + 
	scale_colour_viridis_c() +
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	scale_y_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-14-3.png)<!-- -->

```r
# More accuracy in the middle, definitely!

# reruns model, based around the same variance/covariance structures but with simulated data.
xsim %>% 
	ggplot(aes(y=weighting)) + 
	geom_point(aes(x=sim.val), col="red", size=2, alpha=0.05) +
	geom_point(data= filter(xsim, name == "sim_1"),
			   aes(x=true.val),
			   col="blue", size=2, alpha=0.9) +	
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	labs(x="Heritability estimate", y=expression("Weighting (1/SE"^2*")"))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-14-4.png)<!-- -->

```r
# Seems that the model is in general simulating the true data ok!
# Main problem is that values are simulating beyond the bounds of possible heritability.
```


------------------------------------------------------------------------

### Preamble on further models

Note that the data was previously not complete enough to examine
interactions between explanatory variables properly. For example, some factors such as heritability type and life stage have only a single level represented within some trait types, and some traits involve entirely one level of another factor, e.g., all estimates for nutrient content, immune response, and photochemistry come from single studies, and thus all have only one of the two heritability types within, thus precluding direct comparisons within these trait types.

Thus, to better compare interactions, we conduct separate analyses on subsetted data.

## Which predictor variable first?

To guide our analysis further, we determine which additional factors can be examined in combination with trait type, and compare the sample sizes of each comparison to determine the order of the sub-analyses.


```r
p1 <- data %>% table.plot("stage", "trait", .) # singular for photochemistry, immune response
p_stage <- table.plot("stage", "trait", dat_s)

p2 <- data %>% table.plot("h2", "trait", .) # singular for photochemistry, nutrient content, immune response
p_h2 <- table.plot("h2", "trait", dat_sh)

p3 <- data %>% table.plot("growth.form", "trait", .) # singular for traits: photochemistry, immune response, morphology, immune response, gene expression; singular for columnar and encrusting growth forms
p_growthform <- table.plot("growth.form", "trait", dat_shg)

p4 <- data %>% 	filter(!is.na(temp.diff)) %>%
	table.plot("temp.manip.plus.others", "trait", .) # singular for traits: photochemistry, immune response, morphology, immune response, gene expression; singular for columnar and encrusting growth forms
p_tempmanip <- table.plot("temp.manip.plus.others", "trait", dat_tm)

require(patchwork)
(p1 | p2) / (p3 | p4)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

We see there there are a number of rows/columns with less than 2 levels (i.e. 'singular' terms), which make it impossible to estimate interaction terms for these factor combinations, due to the singularity. Thus, we subset these rows/columns and plot the remaining data to examine the sample sizes and factor levels remaining.


```r
require(patchwork)
(p_stage | p_h2) / (p_growthform | p_tempmanip)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Both life stage and heritability type use similar levels, and thus can be combined for the smaller sample size of heritability type. Life stage thus must be examined first, then heritability type and life stage can be examined together for the smaller data subset. Growth form also uses an even more narrow range of studies, and thus we can examine it after looking at the latter two.

Finally, temperature manipulation results in a few different trait types being removed (e.g., morphology, symbiont community), and the addition of photochemistry, and the data within each trait type are also dissimilar due to some studies being observational and thust being removed. Thus, this analysis must be done on its own, and factors re-added on a case-by-case basis.

## Sub-analysis A: Model selection on dataset allowing a trait type x **life stage** interaction

### Step 0: View 'complete' dataset


```r
dat_s %>% table.plot("stage", "trait", .)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
# exclude photochemistry, immune response
```

Using the smaller data subset, we can now examine possible interactions between trait type and life stage, and include additive terms for other predictors.

### Step 1: Define beyond optimal model

The beyond-optimal model in this case is a full analysis of trait type x life stage:
$$ logit(h^2) \sim N(\mu, \tau^2) \\
\mu = trait \times stage + h^2type + growth\,form +  \epsilon_{i} $$

### Step 2: Select optimal random effects


```r
m.study.nested <- rma.mv(yi=val.log ~ trait * stage + h2 + growth.form,
						 V=sv.log, random = ~1|study/est.id, 
						 data=dat_s, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.spp.nested <- rma.mv(yi=val.log  ~ trait * stage + h2 + growth.form,
						 V=sv.log, random = ~1|species/est.id, 
						 data=dat_s, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.only.est <- rma.mv(yi=val.log ~ trait * stage + h2 + growth.form,
						 V=sv.log, random = ~1|study/est.id, 
						 data=dat_s, method="REML", tdist=T,
						 sigma2 = c(0, NA)) # same as random = ~1|est.id
m.only.study <- rma.mv(yi=val.log ~ trait * stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_s, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|study
m.only.spp <- rma.mv(yi=val.log ~ trait * stage + h2 + growth.form,
					  V=sv.log, random = ~1|species/est.id, 
					  data=dat_s, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|species
m.no.reff <- rma.mv(yi=val.log ~ trait * stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_s, method="REML", tdist=T,
					  sigma2 = c(0, 0)) # no random effects


TableS3 <- myAIC(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff) %>%
	mutate(model = fct_recode(model,
		"Estimate ID nested in study ID" = "m.study.nested",
		"Estimate ID nested in species" = "m.spp.nested",
		"Study ID only" = "m.only.study",
		"Species only" = "m.only.spp",
		"Estimate ID only" = "m.only.est",
		"No random effect" = "m.no.reff"))
write.csv(TableS3, file="Suppl Tables/TableS3.csv", row.names=F)
TableS3
```

<div class="kable-table">

|model                          | df|     AICc|     ΔAICc|
|:------------------------------|--:|--------:|---------:|
|Estimate ID only               | 22| 127.7203|  0.000000|
|Study ID only                  | 22| 130.9134|  3.193072|
|Estimate ID nested in study ID | 22| 133.9795|  6.259191|
|Estimate ID nested in species  | 22| 135.5664|  7.846154|
|No random effect               | 22| 141.7589| 14.038646|
|Species only                   | 22| 148.0740| 20.353683|

</div>

```r
rm(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff)
```

A random effects structure of estimate ID-only is preferred (i.e. to treat estimate IDs as sample sizes of same units). Next preferred is ΔAICc = 1.34.


### Step 3: Select optimal fixed effects


```r
mod_names <- ls(pattern = "m.fixed.")
do.call("rm", as.list(mod_names))

m.fixed.txs.h.g <- rma.mv(yi=val.log ~ trait*stage + h2 + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.t.s.h.g <- rma.mv(yi=val.log ~ trait + stage + h2 + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))


m.fixed.txs.h <- rma.mv(yi=val.log ~ trait*stage + h2,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.txs.g <- rma.mv(yi=val.log ~ trait*stage + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.txs <- rma.mv(yi=val.log ~ trait*stage,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))


m.fixed.t.s.g <- rma.mv(yi=val.log ~ trait + stage + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.t.s.h <- rma.mv(yi=val.log ~ trait + stage + h2,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.t.h.g <- rma.mv(yi=val.log ~ trait + h2 + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.s.h.g <- rma.mv(yi=val.log ~ stage + h2 + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))


m.fixed.t.s <- rma.mv(yi=val.log ~ trait + stage,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.t.h <- rma.mv(yi=val.log ~ trait + h2,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.t.g <- rma.mv(yi=val.log ~ trait + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.s.h <- rma.mv(yi=val.log ~ stage + h2,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.s.g <- rma.mv(yi=val.log ~ stage + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.h.g <- rma.mv(yi=val.log ~ h2 + growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))


m.fixed.t <- rma.mv(yi=val.log ~ trait,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.s <- rma.mv(yi=val.log ~ stage,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.h <- rma.mv(yi=val.log ~ h2,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.g <- rma.mv(yi=val.log ~ growth.form,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))
m.fixed.null <- rma.mv(yi=val.log ~ 1,
					  random = ~1|study/est.id, V=sv.log,
					  data=dat_s, method="ML", tdist=T, sigma2 = c(0, NA))

mod_names <- ls(pattern = "m.fixed.")
TableS4 <- head(do.call("myAIC", as.list(c(mod_names, getmodel=T)))) %>%
	mutate(model = fct_recode(model,
		"trait x life stage + heritability type" = "m.fixed.txs.h",
		"trait x life stage" = "m.fixed.txs",
		"trait x life stage + growth form" = "m.fixed.txs.g",
		"trait x life stage + heritability type + growth form" = "m.fixed.txs.h.g",
		"trait + life stage + growth form" = "m.fixed.t.s.g",
		"trait + heritability type + growth form" = "m.fixed.t.h.g"))
write.csv(TableS4, file="Suppl Tables/TableS4.csv", row.names=F)
TableS4
```

<div class="kable-table">

|model                                                | df|     AICc|     ΔAICc|
|:----------------------------------------------------|--:|--------:|---------:|
|trait x life stage + heritability type               | 18| 80.14774|  0.000000|
|trait x life stage                                   | 17| 82.18116|  2.033424|
|trait x life stage + growth form                     | 21| 88.98944|  8.841700|
|trait x life stage + heritability type + growth form | 22| 90.13248|  9.984738|
|trait + life stage + growth form                     | 12| 91.85181| 11.704067|
|trait + heritability type + growth form              | 11| 92.39415| 12.246413|

</div>

```r
# trait x life stage preferred!

# Double-check with rma fits (knha method is slightly different from rma.mv with t-distributed CIs):
# m.fixed.txs.h <- rma(yi=val.log ~ trait*stage + h2, 
# 					 vi=sv.log, data=dat_s, method="ML", test="knha")
# m.fixed.txs.g <- rma(yi=val.log ~ trait*stage + growth.form, 
# 					 vi=sv.log, data=dat_s, method="ML", test="knha")
# m.fixed.txs <- rma(yi=val.log ~ trait*stage, 
# 					 vi=sv.log, data=dat_s, method="ML", test="knha")
# myAIC(m.fixed.txs.h, m.fixed.txs.g, m.fixed.txs)
# model of trait x life stage + heritability type is still preferred!
```

The trait x life stage + heritability fixed effect structure is optimal. The next closest model drops the effect of heritability type, and is over 2 ΔAICc away (ΔAICc = 2.6).

### Step 4: Fit final model of life stage/heritability type


```r
final.model.A <- rma(yi=val.log ~ trait*stage + h2, vi=sv.log, 
					 data=dat_s, method="REML", test="knha")
final.model.A
```

```

Mixed-Effects Model (k = 74; tau^2 estimator: REML)

tau^2 (estimated amount of residual heterogeneity):     0.0382 (SE = 0.0173)
tau (square root of estimated tau^2 value):             0.1953
I^2 (residual heterogeneity / unaccounted variability): 46.98%
H^2 (unaccounted variability / sampling variability):   1.89
R^2 (amount of heterogeneity accounted for):            78.27%

Test for Residual Heterogeneity:
QE(df = 55) = 96.0189, p-val = 0.0005

Test of Moderators (coefficients 2:19):
F(df1 = 18, df2 = 55) = 7.8780, p-val < .0001

Model Results:

                                       estimate      se     tval    pval 
intrcpt                                 -0.2731  0.5329  -0.5126  0.6103 
traitgrowth                              0.4395  0.1690   2.6012  0.0119 
traitnutrient content                    0.0561  0.5609   0.1001  0.9206 
traitbleaching                          -0.3175  0.5803  -0.5470  0.5866 
traitsymbiont community                  0.1795  0.5609   0.3200  0.7501 
traitmorphology                          0.1718  0.5757   0.2984  0.7665 
traitsurvival                            0.2574  0.5368   0.4794  0.6335 
stagejuvenile                            0.1620  0.1729   0.9374  0.3526 
stageadult                              -0.6607  0.5544  -1.1917  0.2385 
h2narrow                                -0.2716  0.1195  -2.2720  0.0270 
traitgrowth:stagejuvenile               -1.4359  0.5894  -2.4364  0.0181 
traitbleaching:stagejuvenile            -0.6695  0.3533  -1.8951  0.0633 
traitsymbiont community:stagejuvenile   -0.5102  0.6872  -0.7424  0.4610 
traitmorphology:stagejuvenile           -0.4727  0.3466  -1.3640  0.1781 
traitnutrient content:stageadult        -0.1588  0.5932  -0.2678  0.7899 
traitbleaching:stageadult                1.1584  0.6216   1.8636  0.0677 
traitsymbiont community:stageadult       0.6817  0.6128   1.1125  0.2708 
traitmorphology:stageadult               0.3645  0.6503   0.5605  0.5774 
traitsurvival:stageadult                 0.3266  0.5682   0.5748  0.5678 
                                         ci.lb    ci.ub 
intrcpt                                -1.3410   0.7947    
traitgrowth                             0.1009   0.7781  * 
traitnutrient content                  -1.0679   1.1802    
traitbleaching                         -1.4805   0.8456    
traitsymbiont community                -0.9446   1.3036    
traitmorphology                        -0.9819   1.3255    
traitsurvival                          -0.8184   1.3331    
stagejuvenile                          -0.1844   0.5084    
stageadult                             -1.7718   0.4504    
h2narrow                               -0.5112  -0.0320  * 
traitgrowth:stagejuvenile              -2.6170  -0.2548  * 
traitbleaching:stagejuvenile           -1.3776   0.0385  . 
traitsymbiont community:stagejuvenile  -1.8873   0.8670    
traitmorphology:stagejuvenile          -1.1672   0.2218    
traitnutrient content:stageadult       -1.3476   1.0299    
traitbleaching:stageadult              -0.0873   2.4042  . 
traitsymbiont community:stageadult     -0.5463   1.9098    
traitmorphology:stageadult             -0.9388   1.6678    
traitsurvival:stageadult               -0.8121   1.4653    

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
save(final.model.A, file="Models/final.model.A.Rdata")

final.model.A$I
```

```
[1] 46.97796
```

```r
final.model.A$k - final.model.A$p # df for parameter tests
```

```
[1] 55
```

Estimate ID-only, so all residual heterogeneity is: 
$I^2$ = 46.10%, N = 74

### Step 5: Check final model diagnostics

#### Funnel plot


```r
x <- resid(final.model.A)
hist(x) # looks approx. normal, except for one strong outlier at -1.5
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

```r
par(mfrow=c(2,1))
funnel(final.model.A)
update(final.model.A, .~., data=dat_s[x!=0,]) %>%
	funnel(., yaxis="vinv")
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-21-2.png)<!-- -->

```r
par(mfrow=c(1,1))
```

#### Fail-safe number


```r
fsn(yi=val.log, vi=sv.log, data=dat_s, type="Rosenberg") # N=513 to get no significance (but assumes fixed effects)
```

```

Fail-safe N Calculation Using the Rosenberg Approach

Average Effect Size:        -0.0770
Observed Significance Level: <.0001
Target Significance Level:   0.05

Fail-safe N: 505
```

```r
# Number of studies x5 + 10:
length(unique(dat_s$study))*5 + 10
```

```
[1] 100
```


#### Cook's distance


```r
cdist <- cooks.distance(final.model.A, reestimate=FALSE)
plot(cdist)
cdist.i <- which(cdist >2)
points(cdist.i, cdist[cdist.i], col="red", pch=16)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
data[cdist.i,] %>% select(study, species, val, trait, stage, h2) %>%
	cbind(., cooks = round(cdist[cdist.i], 2))
```

<div class="kable-table">

|   |study                |species             |  val|trait              |stage    |h2     | cooks|
|:--|:--------------------|:-------------------|----:|:------------------|:--------|:------|-----:|
|3  |Quigley et al. 2020  |Acropora spathulata | 0.93|survival           |juvenile |narrow |  2.87|
|19 |Wright et al. 2019   |Acropora millepora  | 0.16|immune response    |adult    |broad  |  3.94|
|27 |Manzello et al. 2019 |Orbicella faveolata | 0.73|symbiont community |adult    |broad  |  3.87|

</div>

```r
update(final.model.A, .~., data=dat_s[-cdist.i,]) # effect of juv:growth n.s.
```

```

Mixed-Effects Model (k = 71; tau^2 estimator: REML)

tau^2 (estimated amount of residual heterogeneity):     0.0312 (SE = 0.0167)
tau (square root of estimated tau^2 value):             0.1767
I^2 (residual heterogeneity / unaccounted variability): 37.71%
H^2 (unaccounted variability / sampling variability):   1.61
R^2 (amount of heterogeneity accounted for):            80.78%

Test for Residual Heterogeneity:
QE(df = 52) = 77.6553, p-val = 0.0121

Test of Moderators (coefficients 2:19):
F(df1 = 18, df2 = 52) = 7.7445, p-val < .0001

Model Results:

                                       estimate      se     tval    pval 
intrcpt                                 -0.1933  0.5276  -0.3664  0.7156 
traitgrowth                              0.4393  0.1631   2.6927  0.0095 
traitnutrient content                   -0.0240  0.5527  -0.0433  0.9656 
traitbleaching                          -0.3937  0.5722  -0.6882  0.4944 
traitsymbiont community                  0.1831  0.5498   0.3331  0.7404 
traitmorphology                          0.0920  0.5650   0.1628  0.8713 
traitsurvival                            0.2231  0.5279   0.4227  0.6743 
stagejuvenile                           -0.0201  0.1883  -0.1066  0.9155 
stageadult                              -0.7347  0.5481  -1.3405  0.1859 
h2narrow                                -0.3514  0.1316  -2.6704  0.0101 
traitgrowth:stagejuvenile               -1.2874  0.5805  -2.2179  0.0310 
traitbleaching:stagejuvenile            -0.4081  0.3714  -1.0989  0.2769 
traitsymbiont community:stagejuvenile   -0.3317  0.6814  -0.4868  0.6285 
traitmorphology:stagejuvenile           -0.6400  0.4486  -1.4267  0.1596 
traitnutrient content:stageadult        -0.0860  0.5831  -0.1476  0.8833 
traitbleaching:stageadult                1.0240  0.6382   1.6045  0.1147 
traitsymbiont community:stageadult       0.6723  0.5960   1.1279  0.2645 
traitmorphology:stageadult               0.4698  0.6392   0.7349  0.4657 
traitsurvival:stageadult                 0.3567  0.5575   0.6397  0.5251 
                                         ci.lb    ci.ub 
intrcpt                                -1.2520   0.8654     
traitgrowth                             0.1119   0.7667  ** 
traitnutrient content                  -1.1330   1.0851     
traitbleaching                         -1.5419   0.7544     
traitsymbiont community                -0.9201   1.2864     
traitmorphology                        -1.0417   1.2257     
traitsurvival                          -0.8362   1.2825     
stagejuvenile                          -0.3979   0.3578     
stageadult                             -1.8345   0.3651     
h2narrow                               -0.6155  -0.0874   * 
traitgrowth:stagejuvenile              -2.4522  -0.1226   * 
traitbleaching:stagejuvenile           -1.1533   0.3371     
traitsymbiont community:stagejuvenile  -1.6991   1.0357     
traitmorphology:stagejuvenile          -1.5402   0.2601     
traitnutrient content:stageadult       -1.2561   1.0840     
traitbleaching:stageadult              -0.2566   2.3046     
traitsymbiont community:stageadult     -0.5237   1.8683     
traitmorphology:stageadult             -0.8130   1.7525     
traitsurvival:stageadult               -0.7621   1.4754     

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
update(final.model.A, .~., data=dat_s[-cdist.i[1],]) # immune : temp interaction remains
```

```

Mixed-Effects Model (k = 73; tau^2 estimator: REML)

tau^2 (estimated amount of residual heterogeneity):     0.0367 (SE = 0.0175)
tau (square root of estimated tau^2 value):             0.1916
I^2 (residual heterogeneity / unaccounted variability): 42.32%
H^2 (unaccounted variability / sampling variability):   1.73
R^2 (amount of heterogeneity accounted for):            78.12%

Test for Residual Heterogeneity:
QE(df = 54) = 88.0509, p-val = 0.0023

Test of Moderators (coefficients 2:19):
F(df1 = 18, df2 = 54) = 7.7700, p-val < .0001

Model Results:

                                       estimate      se     tval    pval 
intrcpt                                 -0.1934  0.5267  -0.3673  0.7149 
traitgrowth                              0.4394  0.1659   2.6495  0.0106 
traitnutrient content                   -0.0236  0.5537  -0.0426  0.9662 
traitbleaching                          -0.3965  0.5729  -0.6920  0.4919 
traitsymbiont community                  0.1802  0.5522   0.3264  0.7454 
traitmorphology                          0.0921  0.5678   0.1622  0.8718 
traitsurvival                            0.2155  0.5291   0.4073  0.6854 
stagejuvenile                           -0.0252  0.1976  -0.1273  0.8992 
stageadult                              -0.7392  0.5477  -1.3496  0.1828 
h2narrow                                -0.3513  0.1251  -2.8084  0.0069 
traitgrowth:stagejuvenile               -1.2859  0.5855  -2.1963  0.0324 
traitbleaching:stagejuvenile            -0.4027  0.3763  -1.0702  0.2893 
traitsymbiont community:stagejuvenile   -0.3237  0.6845  -0.4729  0.6382 
traitmorphology:stagejuvenile           -0.2267  0.3643  -0.6224  0.5363 
traitnutrient content:stageadult        -0.0805  0.5852  -0.1376  0.8911 
traitbleaching:stageadult                1.2575  0.6138   2.0488  0.0454 
traitsymbiont community:stageadult       0.6799  0.6024   1.1287  0.2640 
traitmorphology:stageadult               0.4721  0.6421   0.7352  0.4654 
traitsurvival:stageadult                 0.3676  0.5598   0.6567  0.5142 
                                         ci.lb    ci.ub 
intrcpt                                -1.2493   0.8625     
traitgrowth                             0.1069   0.7720   * 
traitnutrient content                  -1.1338   1.0866     
traitbleaching                         -1.5451   0.7522     
traitsymbiont community                -0.9269   1.2873     
traitmorphology                        -1.0464   1.2305     
traitsurvival                          -0.8453   1.2763     
stagejuvenile                          -0.4214   0.3711     
stageadult                             -1.8374   0.3589     
h2narrow                               -0.6021  -0.1005  ** 
traitgrowth:stagejuvenile              -2.4598  -0.1121   * 
traitbleaching:stagejuvenile           -1.1572   0.3518     
traitsymbiont community:stagejuvenile  -1.6961   1.0487     
traitmorphology:stagejuvenile          -0.9571   0.5036     
traitnutrient content:stageadult       -1.2538   1.0928     
traitbleaching:stageadult               0.0269   2.4880   * 
traitsymbiont community:stageadult     -0.5278   1.8875     
traitmorphology:stageadult             -0.8153   1.7595     
traitsurvival:stageadult               -0.7547   1.4899     

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
update(final.model.A, .~., data=dat_s[-cdist.i[2],]) # no interaction after removal!
```

```

Mixed-Effects Model (k = 73; tau^2 estimator: REML)

tau^2 (estimated amount of residual heterogeneity):     0.0385 (SE = 0.0175)
tau (square root of estimated tau^2 value):             0.1963
I^2 (residual heterogeneity / unaccounted variability): 47.24%
H^2 (unaccounted variability / sampling variability):   1.90
R^2 (amount of heterogeneity accounted for):            78.55%

Test for Residual Heterogeneity:
QE(df = 54) = 95.1290, p-val = 0.0005

Test of Moderators (coefficients 2:19):
F(df1 = 18, df2 = 54) = 8.0025, p-val < .0001

Model Results:

                                       estimate      se     tval    pval 
intrcpt                                 -0.2318  0.5313  -0.4363  0.6644 
traitgrowth                              0.4396  0.1683   2.6122  0.0116 
traitnutrient content                    0.0148  0.5593   0.0265  0.9789 
traitbleaching                          -0.3590  0.5786  -0.6205  0.5376 
traitsymbiont community                  0.1793  0.5581   0.3213  0.7492 
traitmorphology                          0.1305  0.5741   0.2272  0.8211 
traitsurvival                            0.2342  0.5345   0.4382  0.6630 
stagejuvenile                            0.1664  0.1725   0.9641  0.3393 
stageadult                              -0.7023  0.5527  -1.2706  0.2093 
h2narrow                                -0.3129  0.1241  -2.5219  0.0147 
traitgrowth:stagejuvenile               -1.4603  0.5869  -2.4881  0.0160 
traitbleaching:stagejuvenile            -0.6325  0.3534  -1.7896  0.0791 
traitsymbiont community:stagejuvenile   -0.5143  0.6838  -0.7521  0.4552 
traitmorphology:stagejuvenile           -0.8264  0.4563  -1.8112  0.0757 
traitnutrient content:stageadult        -0.1171  0.5914  -0.1981  0.8437 
traitbleaching:stageadult                1.2096  0.6202   1.9503  0.0563 
traitsymbiont community:stageadult       0.6822  0.6100   1.1183  0.2684 
traitmorphology:stageadult               0.4207  0.6491   0.6482  0.5196 
traitsurvival:stageadult                 0.3500  0.5657   0.6186  0.5388 
                                         ci.lb    ci.ub 
intrcpt                                -1.2970   0.8334    
traitgrowth                             0.1022   0.7769  * 
traitnutrient content                  -1.1064   1.1361    
traitbleaching                         -1.5189   0.8010    
traitsymbiont community                -0.9397   1.2983    
traitmorphology                        -1.0205   1.2814    
traitsurvival                          -0.8373   1.3057    
stagejuvenile                          -0.1796   0.5123    
stageadult                             -1.8105   0.4059    
h2narrow                               -0.5617  -0.0642  * 
traitgrowth:stagejuvenile              -2.6371  -0.2836  * 
traitbleaching:stagejuvenile           -1.3411   0.0761  . 
traitsymbiont community:stagejuvenile  -1.8852   0.8566    
traitmorphology:stagejuvenile          -1.7413   0.0884  . 
traitnutrient content:stageadult       -1.3028   1.0685    
traitbleaching:stageadult              -0.0338   2.4531  . 
traitsymbiont community:stageadult     -0.5408   1.9053    
traitmorphology:stageadult             -0.8806   1.7220    
traitsurvival:stageadult               -0.7843   1.4842    

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

#### Residuals and simulated data


```r
plot(fitted(final.model.A), rstandard(final.model.A)$z)
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

```r
plot(residuals(final.model.A, type="pearson"))
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-24-2.png)<!-- -->

```r
# Not amazing, not terrible

# Simulate new data given the parameter values of the data and same variance structure, plot on real data to see if it is similar
xsim <- simulate(final.model.A, 100) %>%
	mutate(true.val = dat_s$val.log,
		   weighting = 1/dat_s$sv.log) %>%
	pivot_longer(cols=-c(true.val, weighting)) %>%
	select(name, true.val, sim.val = value, weighting)
xsim %>%
	ggplot(aes(x=true.val, y=sim.val)) +
	geom_point(aes(col=weighting), alpha=0.05) +
	geom_abline(slope=1, intercept=0, linetype="dashed") + 
	scale_colour_viridis_c() +
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	scale_y_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-24-3.png)<!-- -->

```r
# reruns model, based around the same variance/covariance structures but with simulated data.
xsim %>% 
	ggplot(aes(y=weighting)) + 
	geom_point(aes(x=sim.val), col="red", size=2, alpha=0.05) +
	geom_point(data= filter(xsim, name == "sim_1"),
			   aes(x=true.val),
			   col="blue", size=2, alpha=0.9) +	
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	labs(x="Heritability estimate", y=expression("Weighting (1/SE"^2*")"))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-24-4.png)<!-- -->

```r
# Seems that the model is in general simulating the true data ok!
# Main problem is that values are simulating beyond the bounds of possible heritability.
```

------------------------------------------------------------------------

## Sub-analysis B: Model selection on dataset allowing a trait type x **heritability type** interaction

Main result: trait type x heritability type interaction supported. Most
of the interaction is likely driven by the gene expression trait level,
where one estimate of narrow-sense h2 is higher than all broad-sense
estimates of H2.

### Step 0: View 'complete' dataset


```r
dat_sh %>% table.plot("h2", "trait", .) # subset allows examination of h2 x trait type
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

```r
dat_sh %>% table.plot("stage", "trait", .) # can also look at life stage interaction too!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-25-2.png)<!-- -->

```r
dat_sh %>% table.plot("stage", "h2", .)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-25-3.png)<!-- -->

```r
dat_sh %>% table.plot("stage", "growth.form", .) # some bad levels for growth form!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-25-4.png)<!-- -->

The data subset allows for the examination of h2 x trait type interaction, as well as the trait x life stage interaction seen previously!


### Step 1: Define beyond optimal model

The beyond-optimal model in this case is a full analysis of trait type x heritability type + trait type x life stage, plus additive effects:
$$ logit(h^2) \sim N(\mu, \tau^2) \\
\mu = trait \times stage + trait \times h^2type + growth\,form +  \epsilon_{i} $$

### Step 2: Select optimal random effects


```r
m.study.nested <- rma.mv(yi=val.log ~ trait * h2 + trait * stage + growth.form,
						 V=sv.log, random = ~1|study/est.id, 
						 data=dat_sh, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.spp.nested <- rma.mv(yi=val.log ~ trait * h2 + trait * stage + growth.form,
						 V=sv.log, random = ~1|species/est.id, 
						 data=dat_sh, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.only.est <- rma.mv(yi=val.log ~ trait * h2 + trait * stage + growth.form,
						 V=sv.log, random = ~1|study/est.id, 
						 data=dat_sh, method="REML", tdist=T,
						 sigma2 = c(0, NA)) # same as random = ~1|est.id
m.only.study <- rma.mv(yi=val.log ~ trait * h2 + trait * stage + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_sh, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|study
m.only.spp <- rma.mv(yi=val.log ~ trait * h2 + trait * stage + growth.form,
					  V=sv.log, random = ~1|species/est.id, 
					  data=dat_sh, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|species
m.no.reff <- rma.mv(yi=val.log ~ trait * h2 + trait * stage + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_sh, method="REML", tdist=T,
					  sigma2 = c(0, 0)) # no random effects


TableS5 <- myAIC(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff) %>%
	mutate(model = fct_recode(model,
		"Estimate ID nested in study ID" = "m.study.nested",
		"Estimate ID nested in species" = "m.spp.nested",
		"Study ID only" = "m.only.study",
		"Species only" = "m.only.spp",
		"Estimate ID only" = "m.only.est",
		"No random effect" = "m.no.reff"))
write.csv(TableS5, file="Suppl Tables/TableS5.csv", row.names=F)
TableS5
```

<div class="kable-table">

|model                          | df|     AICc|      ΔAICc|
|:------------------------------|--:|--------:|----------:|
|Study ID only                  | 23| 161.1380|  0.0000000|
|Estimate ID only               | 23| 161.1796|  0.0416663|
|No random effect               | 23| 168.1204|  6.9824555|
|Estimate ID nested in study ID | 23| 172.7948| 11.6568354|
|Estimate ID nested in species  | 23| 174.4590| 13.3210780|
|Species only                   | 23| 178.7043| 17.5662968|

</div>

```r
rm(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff)
```
The optimal model is using study ID-only, however the second-most preferred model uses estimate ID only, and ΔAICc = 0.27, so close call either way.
Let's continue with study ID.


### Step 3: Select optimal fixed effects


```r
mod_names <- ls(pattern = "m.fixed.")
do.call("rm", as.list(mod_names))

m.fixed.txh.txs.g <- rma.mv(yi=val.log ~ trait*h2 + trait*stage + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txh.txs <- rma.mv(yi=val.log ~ trait*h2 + trait*stage,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))

m.fixed.txh.s.g <- rma.mv(yi=val.log ~ trait*h2 + stage + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txh.s <- rma.mv(yi=val.log ~ trait*h2 + stage,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txh.g <- rma.mv(yi=val.log ~ trait*h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txh <- rma.mv(yi=val.log ~ trait*h2,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))

m.fixed.txs.h.g <- rma.mv(yi=val.log ~ trait*stage + h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txs.h <- rma.mv(yi=val.log ~ trait*stage + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txs.g <- rma.mv(yi=val.log ~ trait*stage + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txs <- rma.mv(yi=val.log ~ trait*stage,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))

m.fixed.t.s.h.g <- rma.mv(yi=val.log ~ trait + stage + h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.s.h <- rma.mv(yi=val.log ~ trait + stage + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.s.g <- rma.mv(yi=val.log ~ trait + stage + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.h.g <- rma.mv(yi=val.log ~ trait + h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.s.h.g <- rma.mv(yi=val.log ~ stage + h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))

m.fixed.t.s <- rma.mv(yi=val.log ~ trait + stage,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.h <- rma.mv(yi=val.log ~ trait + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.g <- rma.mv(yi=val.log ~ trait + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.s.h <- rma.mv(yi=val.log ~ stage + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.s.g <- rma.mv(yi=val.log ~ stage + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.h.g <- rma.mv(yi=val.log ~ h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))


m.fixed.t <- rma.mv(yi=val.log ~ trait,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.s <- rma.mv(yi=val.log ~ stage,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.h <- rma.mv(yi=val.log ~ h2,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.g <- rma.mv(yi=val.log ~ growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.null <- rma.mv(yi=val.log ~ 1,
							V=sv.log, random = ~1|study/est.id, data=dat_sh,
							method="ML", tdist=T, sigma2 = c(NA, 0))


mod_names <- ls(pattern = "m.fixed.")
TableS6 <- head(do.call("myAIC", as.list(c(mod_names, getmodel=T)))) %>%
	mutate(model = fct_recode(model,
		"trait x life stage" = "m.fixed.txs",
		"trait x life stage + heritability type" = "m.fixed.txs.h",
		"trait x life stage + trait x heritability type" = "m.fixed.txh.txs",
		"trait x heritability type + life stage" = "m.fixed.txh.s",
		"trait + heritability type" = "m.fixed.txh",
		"trait x life stage + growth form" = "m.fixed.txs.g"))
write.csv(TableS6, file="Suppl Tables/TableS6.csv", row.names=F)
TableS6
```

<div class="kable-table">

|model                                          | df|     AICc|     ΔAICc|
|:----------------------------------------------|--:|--------:|---------:|
|trait x life stage                             | 15| 78.15443|  0.000000|
|trait x life stage + heritability type         | 16| 80.83149|  2.677067|
|trait x life stage + trait x heritability type | 19| 84.42971|  6.275287|
|trait + heritability type                      | 11| 87.69613|  9.541702|
|trait x heritability type + life stage         | 13| 89.18243| 11.028005|
|trait x life stage + growth form               | 19| 91.08915| 12.934727|

</div>

A model of trait type x life stage only is preferred! The next closest model is
one with trait type x life stage + heritability type (ΔAICc = 2.15).

### Step 4: Fit final model of heritability type


```r
final.model.B <- rma.mv(yi=val.log ~ trait*stage, 
					  V=sv.log, random = ~1|study/est.id,
					  data=dat_sh, method="REML", tdist=T, sigma2 = c(NA, 0))
rma(yi=val.log ~ trait*stage, 
					  vi=sv.log, random = ~1|study,
					  data=dat_sh, method="REML", test="knha")
```

```

Mixed-Effects Model (k = 67; tau^2 estimator: REML)

tau^2 (estimated amount of residual heterogeneity):     0.0506 (SE = 0.0214)
tau (square root of estimated tau^2 value):             0.2249
I^2 (residual heterogeneity / unaccounted variability): 59.08%
H^2 (unaccounted variability / sampling variability):   2.44
R^2 (amount of heterogeneity accounted for):            71.18%

Test for Residual Heterogeneity:
QE(df = 51) = 98.0787, p-val < .0001

Test of Moderators (coefficients 2:16):
F(df1 = 15, df2 = 51) = 6.3854, p-val < .0001

Model Results:

                                       estimate      se     tval    pval 
intrcpt                                 -0.5447  0.5338  -1.0205  0.3123 
traitgrowth                              0.4409  0.1788   2.4665  0.0170 
traitbleaching                          -0.0513  0.5871  -0.0874  0.9307 
traitsymbiont community                  0.1741  0.5801   0.3001  0.7653 
traitmorphology                          0.4434  0.5869   0.7554  0.4535 
traitsurvival                            0.3882  0.5499   0.7060  0.4834 
stagejuvenile                            0.1248  0.1932   0.6457  0.5214 
stageadult                              -0.3981  0.5575  -0.7140  0.4785 
traitgrowth:stagejuvenile               -1.2857  0.6089  -2.1115  0.0396 
traitbleaching:stagejuvenile            -0.9029  0.3637  -2.4826  0.0164 
traitsymbiont community:stagejuvenile   -0.4675  0.7096  -0.6588  0.5130 
traitmorphology:stagejuvenile           -0.6348  0.3765  -1.6864  0.0978 
traitbleaching:stageadult                0.8123  0.6270   1.2954  0.2010 
traitsymbiont community:stageadult       0.6961  0.6417   1.0848  0.2831 
traitmorphology:stageadult              -0.0036  0.6569  -0.0056  0.9956 
traitsurvival:stageadult                 0.2024  0.5845   0.3463  0.7305 
                                         ci.lb    ci.ub 
intrcpt                                -1.6164   0.5269    
traitgrowth                             0.0820   0.7998  * 
traitbleaching                         -1.2299   1.1273    
traitsymbiont community                -0.9906   1.3387    
traitmorphology                        -0.7350   1.6217    
traitsurvival                          -0.7158   1.4923    
stagejuvenile                          -0.2632   0.5127    
stageadult                             -1.5173   0.7212    
traitgrowth:stagejuvenile              -2.5081  -0.0633  * 
traitbleaching:stagejuvenile           -1.6330  -0.1727  * 
traitsymbiont community:stagejuvenile  -1.8920   0.9570    
traitmorphology:stagejuvenile          -1.3906   0.1209  . 
traitbleaching:stageadult              -0.4466   2.0711    
traitsymbiont community:stageadult     -0.5922   1.9844    
traitmorphology:stageadult             -1.3224   1.3151    
traitsurvival:stageadult               -0.9711   1.3760    

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# final.model.B
save(final.model.B, file="Models/final.model.B.Rdata")

mlm.variance.distribution(final.model.B)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-28-1.png)<!-- --><div class="kable-table">

|        | % of total variance|I2    |
|:-------|-------------------:|:-----|
|Level 1 |            25.04703|---   |
|Level 2 |             0.00000|0     |
|Level 3 |            74.95297|74.95 |

</div>

```r
final.model.B$k - final.model.B$p
```

```
[1] 51
```

Without the effect of heritability type, juvenile bleaching also is a significant interaction term!

### Step 5: Check final model diagnostics

#### Funnel plot


```r
x <- resid(final.model.B)
hist(x) # looks approx. normal, except for one strong outlier at -1.5
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

```r
par(mfrow=c(2,1))
funnel(final.model.B)
funnel(final.model.B, yaxis = "vinv")
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-29-2.png)<!-- -->

```r
par(mfrow=c(1,1))

dat_sh[x < -1,] %>% select(study, species, trait, stage, h2, val)
```

<div class="kable-table">

|   |study               |species             |trait    |stage    |h2     |     val|
|:--|:-------------------|:-------------------|:--------|:--------|:------|-------:|
|15 |Quigley et al. 2020 |Acropora spathulata |survival |juvenile |narrow | 0.00038|

</div>

```r
# update(final.model.B, .~., data=dat_sh[x > -1,]) %>%
	# funnel() # Funnel looks good once that one point is removed!
```

#### Fail-safe number and Cook's distance


```r
fsn(yi=val.log, vi=sv.log, data=dat_sh, type="Rosenberg") # N=513 to get no significance (but assumes fixed effects)
```

```

Fail-safe N Calculation Using the Rosenberg Approach

Average Effect Size:        -0.0477
Observed Significance Level: 0.0009
Target Significance Level:   0.05

Fail-safe N: 127
```

```r
# Number of studies x5 + 10:
length(unique(dat_sh$study))*5 + 10
```

```
[1] 100
```

```r
cdist <- cooks.distance(final.model.B, reestimate=FALSE)
plot(cdist)
cdist.i <- which(cdist >2)
points(cdist.i, cdist[cdist.i], col="red", pch=16)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

```r
dat_sh[cdist.i,] %>% select(study, trait, stage, h2, val) %>% mutate(cdist = cdist[cdist.i])
```

<div class="kable-table">

|study                |trait     |stage    |h2     |       val|     cdist|
|:--------------------|:---------|:--------|:------|---------:|---------:|
|Kenkel et al. 2015   |survival  |juvenile |broad  | 0.9400000|  8.992638|
|Quigley et al. 2020  |survival  |juvenile |narrow | 0.9300000| 13.235339|
|Manzello et al. 2019 |bleaching |adult    |broad  | 0.9300000| 60.246805|
|Zhang et al. 2019    |bleaching |larvae   |broad  | 0.4511525|  2.069310|
|Quigley et al. 2020  |bleaching |juvenile |narrow | 0.1500000|  2.488141|
|Császár et al. 2010  |growth    |adult    |broad  | 0.5900000|  6.604630|

</div>

```r
update(final.model.B, .~., data=dat_sh[-cdist.i[1],])
```

```

Multivariate Meta-Analysis Model (k = 66; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0726  0.2694     18     no         study 
sigma^2.2  0.0000  0.0000     66    yes  study/est.id 

Test for Residual Heterogeneity:
QE(df = 50) = 97.6866, p-val < .0001

Test of Moderators (coefficients 2:16):
F(df1 = 15, df2 = 50) = 13.7493, p-val < .0001

Model Results:

                                       estimate      se     tval    pval 
intrcpt                                 -0.4900  0.5603  -0.8744  0.3861 
traitgrowth                              0.5823  0.1601   3.6365  0.0007 
traitbleaching                           0.0216  0.6175   0.0350  0.9722 
traitsymbiont community                  0.1121  0.6178   0.1814  0.8568 
traitmorphology                          0.4845  0.5986   0.8093  0.4222 
traitsurvival                            0.3098  0.5696   0.5439  0.5890 
stagejuvenile                            0.5172  0.2683   1.9274  0.0596 
stageadult                              -0.5949  0.5961  -0.9980  0.3231 
traitgrowth:stagejuvenile               -1.8358  0.6121  -2.9992  0.0042 
traitbleaching:stagejuvenile            -1.1454  0.3094  -3.7023  0.0005 
traitsymbiont community:stagejuvenile   -0.9659  0.7138  -1.3533  0.1821 
traitmorphology:stagejuvenile           -1.0902  0.4303  -2.5337  0.0145 
traitbleaching:stageadult                0.8833  0.6623   1.3336  0.1884 
traitsymbiont community:stageadult       0.6005  0.6652   0.9027  0.3710 
traitmorphology:stageadult               0.1383  0.7203   0.1920  0.8485 
traitsurvival:stageadult                 0.7554  0.6052   1.2482  0.2178 
                                         ci.lb    ci.ub 
intrcpt                                -1.6154   0.6355      
traitgrowth                             0.2607   0.9039  *** 
traitbleaching                         -1.2187   1.2619      
traitsymbiont community                -1.1288   1.3529      
traitmorphology                        -0.7179   1.6869      
traitsurvival                          -0.8343   1.4539      
stagejuvenile                          -0.0218   1.0561    . 
stageadult                             -1.7923   0.6025      
traitgrowth:stagejuvenile              -3.0653  -0.6064   ** 
traitbleaching:stagejuvenile           -1.7668  -0.5240  *** 
traitsymbiont community:stagejuvenile  -2.3996   0.4677      
traitmorphology:stagejuvenile          -1.9545  -0.2260    * 
traitbleaching:stageadult              -0.4470   2.2135      
traitsymbiont community:stageadult     -0.7356   1.9366      
traitmorphology:stageadult             -1.3085   1.5850      
traitsurvival:stageadult               -0.4602   1.9710      

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# dat_sh %>% table.plot("stage", "trait", .) 
```

Plot the coefficients with different model datasets:

```r
dat_sh.out1 <- update(final.model.B, .~., data=dat_sh[-cdist.i[1],])
dat_sh.out2 <- update(final.model.B, .~., data=dat_sh[-cdist.i[2],])
dat_sh.out3 <- update(final.model.B, .~., data=dat_sh[-cdist.i[3],])
dat_sh.out4 <- update(final.model.B, .~., data=dat_sh[-cdist.i[4],])
dat_sh.out5 <- update(final.model.B, .~., data=dat_sh[-cdist.i[5],])
dat_sh.out.all <- update(final.model.B, .~., data=dat_sh[-cdist.i[c(1,2,3,5)],])


rbind(
	data.frame(
		Model = rep(" Original model fit"),
		Parameter = names(coef(final.model.B)), 
		Estimate = coef(final.model.B),
		SE = final.model.B$se),
	data.frame(
		Model = rep(paste0(" without Kenkel et al. 2015 - juvenile survival\nCook's dist = ",round(cdist[cdist.i[1]],1))),
		Parameter = names(coef(dat_sh.out1)), 
		Estimate = coef(dat_sh.out1),
		SE = dat_sh.out1$se),
	data.frame(
		Model = rep(paste0(" without Quigley et al. 2020 - juvenile survival\nCook's dist = ",round(cdist[cdist.i[2]],1))),
		Parameter = names(coef(dat_sh.out2)), 
		Estimate = coef(dat_sh.out2),
		SE = dat_sh.out2$se),
	data.frame(
		Model = rep(paste0(" without Manzello et al. 2019 - adult bleaching\nCook's dist = ",round(cdist[cdist.i[3]],1))),
		Parameter = names(coef(dat_sh.out3)), 
		Estimate = coef(dat_sh.out3),
		SE = dat_sh.out3$se),
	data.frame(
		Model = rep(paste0(" without Quigley et al. 2020 - juvenile bleaching\nCook's dist = ",round(cdist[cdist.i[4]],1))),
		Parameter = names(coef(dat_sh.out4)), 
		Estimate = coef(dat_sh.out4),
		SE = dat_sh.out4$se),
	data.frame(
		Model = rep(paste0(" without Császár et al. 2010 - adult growth\nCook's dist = ",round(cdist[cdist.i[5]],1))),
		Parameter = names(coef(dat_sh.out5)), 
		Estimate = coef(dat_sh.out5),
		SE = dat_sh.out5$se),
	data.frame(
		Model = rep(paste0("without estimates with Cook's distance > 5")),
		Parameter = names(coef(dat_sh.out.all)), 
		Estimate = coef(dat_sh.out.all),
		SE = dat_sh.out.all$se)) %>%
	mutate(Parameter = factor(Parameter),
		   Parameter = fct_relevel(Parameter, "intrcpt")) %>%
	ggplot(aes(x=Estimate, y=Parameter, color=Model)) +
	geom_vline(xintercept=0, linetype="dashed") +
	geom_pointrange(aes(xmin = Estimate - SE, xmax = Estimate + SE), position = position_dodge(width=0.3)) + 
	labs(x="Parameter estimate (log[X+0.2] scale) ± SE") +
	theme(legend.position="top", legend.direction="vertical")
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-31-1.png)<!-- -->


#### Residuals and simulated data


```r
plot(fitted(final.model.B), rstandard(final.model.B)$z)
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

```r
plot(residuals(final.model.B, type="pearson"))
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-32-2.png)<!-- -->

```r
# Not amazing, not terrible

# Simulate new data given the parameter values of the data and same variance structure, plot on real data to see if it is similar
xsim <- simulate(final.model.B, 100) %>%
	mutate(true.val = dat_sh$val.log,
		   weighting = 1/dat_sh$sv.log) %>%
	pivot_longer(cols=-c(true.val, weighting)) %>%
	select(name, true.val, sim.val = value, weighting)
xsim %>%
	ggplot(aes(x=true.val, y=sim.val)) +
	geom_point(aes(col=weighting), alpha=0.05) +
	geom_abline(slope=1, intercept=0, linetype="dashed") + 
	scale_colour_viridis_c() +
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	scale_y_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-32-3.png)<!-- -->

```r
# reruns model, based around the same variance/covariance structures but with simulated data.
xsim %>% 
	ggplot(aes(y=weighting)) + 
	geom_point(aes(x=sim.val), col="red", size=2, alpha=0.05) +
	geom_point(data= filter(xsim, name == "sim_1"),
			   aes(x=true.val),
			   col="blue", size=2, alpha=0.9) +	
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	labs(x="Heritability estimate", y=expression("Weighting (1/SE"^2*")"))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-32-4.png)<!-- -->

```r
# Seems that the model is in general simulating the true data ok!
# Main problem is that values are simulating beyond the bounds of possible heritability.
```

------------------------------------------------------------------------

## Part C: Sub-analysis of trait type x **growth form** interaction

### Step 0: View 'complete' dataset

For this analysis, note that there are only three estimates for
encrusting species (*Montipora capitata*, *Montipora flabellate*,
*Montipora patula*) and only a single estimate for columnar coral
(*Porites evermanni*). Thus, we initially lumped these species into an 'Other' category, but quickly found that they lacked most of the trait and life stage combinations in order to be included in any analysis.


```r
dat_shg %>% table.plot("growth.form", "trait", .)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

```r
dat_shg %>% table.plot("stage", "trait", .) # can also look at life stage interaction too
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-33-2.png)<!-- -->

```r
dat_shg %>% table.plot("stage", "h2", .)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-33-3.png)<!-- -->

```r
dat_shg %>% table.plot("h2", "growth.form", .)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-33-4.png)<!-- -->

### Step 1: Define beyond optimal model

The beyond-optimal model in this case is a full analysis of trait type x growth form + trait type x life stage, plus additive effects:
$$ logit(h^2) \sim N(\mu, \tau^2) \\
\mu = trait \times growth\,form + trait \times stage + h^2type + \epsilon_{i} $$

### Step 2: Select optimal random effects


```r
# growth form without h2:
m.study.nested <- rma.mv(yi=val.log ~ trait * growth.form + trait * stage + h2,
						 V=sv.log, random = ~1|study/est.id,
						 data=dat_shg, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.spp.nested <- rma.mv(yi=val.log  ~ trait * growth.form + trait * stage + h2,
						 V=sv.log, random = ~1|species/est.id,
						 data=dat_shg, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.only.est <- rma.mv(yi=val.log ~ trait * growth.form + trait * stage + h2,
						 V=sv.log, random = ~1|study/est.id,
						 data=dat_shg, method="REML", tdist=T,
						 sigma2 = c(0, NA)) # same as random = ~1|est.id
m.only.study <- rma.mv(yi=val.log ~ trait * growth.form + trait * stage + h2,
					  V=sv.log, random = ~1|study/est.id,
					  data=dat_shg, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|study
m.only.spp <- rma.mv(yi=val.log ~ trait * growth.form + trait * stage + h2,
					  V=sv.log, random = ~1|species/est.id,
					  data=dat_shg, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|species
m.no.reff <- rma.mv(yi=val.log ~ trait * growth.form + trait * stage + h2,
					  V=sv.log, random = ~1|study/est.id,
					  data=dat_shg, method="REML", tdist=T,
					  sigma2 = c(0, 0)) # no random effects


TableS7 <- myAIC(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff) %>%
	mutate(model = fct_recode(model,
		"Estimate ID nested in study ID" = "m.study.nested",
		"Estimate ID nested in species" = "m.spp.nested",
		"Study ID only" = "m.only.study",
		"Species only" = "m.only.spp",
		"Estimate ID only" = "m.only.est",
		"No random effect" = "m.no.reff"))
write.csv(TableS7, file="Suppl Tables/TableS7.csv", row.names=F)
TableS7
```

<div class="kable-table">

|model                          | df|     AICc|      ΔAICc|
|:------------------------------|--:|--------:|----------:|
|Study ID only                  | 18| 120.6457|  0.0000000|
|Estimate ID only               | 18| 121.1527|  0.5070027|
|No random effect               | 18| 122.5166|  1.8708784|
|Estimate ID nested in study ID | 18| 132.7621| 12.1163709|
|Species only                   | 18| 133.0640| 12.4182707|
|Estimate ID nested in species  | 18| 134.2297| 13.5839257|

</div>

```r
rm(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff)
```

Study ID only is preferred, estimate ID is close to optimal as well (ΔAIC=0.5), and null effects structure may also be preferred (ΔAIC=1.9). Continuing with study ID as the random effect...

### Step 3: Select optimal fixed effects



```r
mod_names <- ls(pattern = "m.fixed.")
do.call("rm", as.list(mod_names))

m.fixed.txg.txs.h <- rma.mv(yi=val.log ~ trait*growth.form + trait*stage + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txg.txs <- rma.mv(yi=val.log ~ trait*growth.form + trait*stage,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))

m.fixed.txg.s.h <- rma.mv(yi=val.log ~ trait*growth.form + stage + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txg.s <- rma.mv(yi=val.log ~ trait*growth.form + stage,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txg.h <- rma.mv(yi=val.log ~ trait*growth.form + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txg <- rma.mv(yi=val.log ~ trait*growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))

m.fixed.txs.h.g <- rma.mv(yi=val.log ~ trait*stage + h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txs.h <- rma.mv(yi=val.log ~ trait*stage + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txs.g <- rma.mv(yi=val.log ~ trait*stage + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.txs <- rma.mv(yi=val.log ~ trait*stage,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))

m.fixed.t.s.h.g <- rma.mv(yi=val.log ~ trait + stage + h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.s.h <- rma.mv(yi=val.log ~ trait + stage + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.s.g <- rma.mv(yi=val.log ~ trait + stage + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.h.g <- rma.mv(yi=val.log ~ trait + h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.s.h.g <- rma.mv(yi=val.log ~ stage + h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))

m.fixed.t.s <- rma.mv(yi=val.log ~ trait + stage,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.h <- rma.mv(yi=val.log ~ trait + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.t.g <- rma.mv(yi=val.log ~ trait + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.s.h <- rma.mv(yi=val.log ~ stage + h2,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.s.g <- rma.mv(yi=val.log ~ stage + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.h.g <- rma.mv(yi=val.log ~ h2 + growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))


m.fixed.t <- rma.mv(yi=val.log ~ trait,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.s <- rma.mv(yi=val.log ~ stage,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.h <- rma.mv(yi=val.log ~ h2,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.g <- rma.mv(yi=val.log ~ growth.form,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))
m.fixed.null <- rma.mv(yi=val.log ~ 1,
							V=sv.log, random = ~1|study/est.id, data=dat_shg,
							method="ML", tdist=T, sigma2 = c(NA, 0))


mod_names <- ls(pattern = "m.fixed.")
TableS8 <- head(do.call("myAIC", as.list(c(mod_names, getmodel=T)))) %>%
	mutate(model = fct_recode(model,
		"trait x life stage" = "m.fixed.txs",
		"trait x life stage + heritability type" = "m.fixed.txs.h",
		"trait x life stage + trait x growth form + heritability type" = "m.fixed.txg.txs.h",
		"trait x life stage + trait x growth form" = "m.fixed.txg.txs",
		"trait x life stage + growth form" = "m.fixed.txs.g",
		"trait x life stage + heritability type + growth form" = "m.fixed.txs.h.g"))
write.csv(TableS8, file="Suppl Tables/TableS8.csv", row.names=F)
TableS8
```

<div class="kable-table">

|model                                                        | df|     AICc|    ΔAICc|
|:------------------------------------------------------------|--:|--------:|--------:|
|trait x life stage                                           | 12| 57.98042| 0.000000|
|trait x life stage + heritability type                       | 13| 60.54629| 2.565869|
|trait x life stage + growth form                             | 14| 64.92559| 6.945171|
|trait x life stage + trait x growth form + heritability type | 18| 64.99195| 7.011524|
|trait x life stage + trait x growth form                     | 17| 65.29626| 7.315839|
|trait x life stage + heritability type + growth form         | 15| 67.87557| 9.895152|

</div>

### Step 4: Fit final model


```r
final.model.C <- rma.mv(yi=val.log ~ trait*stage, 
					  V=sv.log, random = ~1|study/est.id,
					  data=dat_shg, method="REML", tdist=T, sigma2 = c(NA, 0))
final.model.C
```

```

Multivariate Meta-Analysis Model (k = 54; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0529  0.2299     17     no         study 
sigma^2.2  0.0000  0.0000     54    yes  study/est.id 

Test for Residual Heterogeneity:
QE(df = 41) = 93.9742, p-val < .0001

Test of Moderators (coefficients 2:13):
F(df1 = 12, df2 = 41) = 23.6844, p-val < .0001

Model Results:

                                       estimate      se     tval    pval 
intrcpt                                 -0.6172  0.1912  -3.2273  0.0025 
traitnutrient content                    0.4909  0.2465   1.9919  0.0531 
traitbleaching                           0.1467  0.2912   0.5039  0.6170 
traitsymbiont community                  0.2457  0.3067   0.8011  0.4277 
traitsurvival                            0.4578  0.1300   3.5205  0.0011 
stagejuvenile                           -0.5560  0.2683  -2.0725  0.0445 
stageadult                               0.1402  0.1997   0.7022  0.4865 
traitbleaching:stagejuvenile            -0.2689  0.3408  -0.7891  0.4346 
traitsymbiont community:stagejuvenile    0.1181  0.7130   0.1656  0.8693 
traitsurvival:stagejuvenile              0.8090  0.1631   4.9600  <.0001 
traitnutrient content:stageadult        -0.7353  0.2389  -3.0780  0.0037 
traitbleaching:stageadult                0.1797  0.3139   0.5725  0.5701 
traitsymbiont community:stageadult      -0.1109  0.3336  -0.3325  0.7412 
                                         ci.lb    ci.ub 
intrcpt                                -1.0034  -0.2310   ** 
traitnutrient content                  -0.0068   0.9887    . 
traitbleaching                         -0.4413   0.7347      
traitsymbiont community                -0.3736   0.8650      
traitsurvival                           0.1952   0.7204   ** 
stagejuvenile                          -1.0977  -0.0142    * 
stageadult                             -0.2630   0.5434      
traitbleaching:stagejuvenile           -0.9572   0.4193      
traitsymbiont community:stagejuvenile  -1.3218   1.5580      
traitsurvival:stagejuvenile             0.4796   1.1384  *** 
traitnutrient content:stageadult       -1.2177  -0.2528   ** 
traitbleaching:stageadult              -0.4542   0.8136      
traitsymbiont community:stageadult     -0.7847   0.5628      

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
save(final.model.C, file="Models/final.model.C.Rdata")

mlm.variance.distribution(final.model.C)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-36-1.png)<!-- --><div class="kable-table">

|        | % of total variance|I2    |
|:-------|-------------------:|:-----|
|Level 1 |            20.66957|---   |
|Level 2 |             0.00000|0     |
|Level 3 |            79.33043|79.33 |

</div>

```r
final.model.C$k - final.model.C$p
```

```
[1] 41
```


A final fixed effect structure of ~trait type x life stage + heritability type (for estimate-ID only random effect) or ~trait type x life stage (for study-ID only random effect).


### Step 5: Check final model diagnostics

#### Funnel plot


```r
x <- resid(final.model.C)
hist(x) # looks approx. normal, except for one strong outlier at -1.5
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-37-1.png)<!-- -->

```r
par(mfrow=c(2,1))
funnel(final.model.C)
funnel(final.model.C, yaxis = "vinv")
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-37-2.png)<!-- -->

```r
par(mfrow=c(1,1))

dat_shg[x < -1,] %>% select(study, species, trait, stage, h2, val)
```

<div class="kable-table">

|   |study               |species             |trait    |stage    |h2     |     val|
|:--|:-------------------|:-------------------|:--------|:--------|:------|-------:|
|15 |Quigley et al. 2020 |Acropora spathulata |survival |juvenile |narrow | 0.00038|

</div>

```r
update(final.model.C, .~., data=dat_shg[x > -1,]) %>%
	funnel()
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-37-3.png)<!-- -->

Funnel looks good once that one category with n=1 is removed!

#### Fail-safe number  Cook's distance


```r
fsn(yi=val.log, vi=sv.log, data=dat_shg, type="Rosenberg") # N=513 to get no significance (but assumes fixed effects)
```

```

Fail-safe N Calculation Using the Rosenberg Approach

Average Effect Size:        -0.0512
Observed Significance Level: 0.0004
Target Significance Level:   0.05

Fail-safe N: 122
```

```r
# Number of studies x5 + 10:
length(unique(dat_shg$study))*5 + 10
```

```
[1] 95
```

```r
cdist <- cooks.distance(final.model.C, reestimate=FALSE)
plot(cdist)
cdist.i <- which(cdist >2)
points(cdist.i, cdist[cdist.i], col="red", pch=16)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-38-1.png)<!-- -->

```r
dat_shg[cdist.i,] %>% select(study, trait, stage, h2, val) %>% mutate(cdist = cdist[cdist.i])
```

<div class="kable-table">

|study                |trait     |stage    |h2     |       val|     cdist|
|:--------------------|:---------|:--------|:------|---------:|---------:|
|Kenkel et al. 2015   |survival  |juvenile |broad  | 0.9400000|  9.213022|
|Quigley et al. 2020  |survival  |juvenile |narrow | 0.9300000| 13.513205|
|Manzello et al. 2019 |bleaching |adult    |broad  | 0.9300000| 58.185574|
|Zhang et al. 2019    |bleaching |larvae   |broad  | 0.4511525|  2.085304|
|Quigley et al. 2020  |bleaching |juvenile |narrow | 0.1500000|  2.487176|

</div>

```r
update(final.model.C, .~., data=dat_shg[-cdist.i[3],])
```

```

Multivariate Meta-Analysis Model (k = 53; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0393  0.1983     17     no         study 
sigma^2.2  0.0000  0.0000     53    yes  study/est.id 

Test for Residual Heterogeneity:
QE(df = 40) = 80.6448, p-val = 0.0001

Test of Moderators (coefficients 2:13):
F(df1 = 12, df2 = 40) = 23.0564, p-val < .0001

Model Results:

                                       estimate      se     tval    pval 
intrcpt                                 -0.5565  0.1800  -3.0911  0.0036 
traitnutrient content                    0.4251  0.2385   1.7826  0.0822 
traitbleaching                           0.0809  0.2844   0.2844  0.7776 
traitsymbiont community                  0.1907  0.2877   0.6629  0.5112 
traitsurvival                            0.4162  0.1280   3.2522  0.0023 
stagejuvenile                           -0.6165  0.2471  -2.4954  0.0168 
stageadult                               0.0554  0.1862   0.2974  0.7677 
traitbleaching:stagejuvenile            -0.2038  0.3350  -0.6084  0.5464 
traitsymbiont community:stagejuvenile    0.1884  0.7046   0.2674  0.7906 
traitsurvival:stagejuvenile              0.8508  0.1615   5.2696  <.0001 
traitnutrient content:stageadult        -0.7129  0.2317  -3.0766  0.0038 
traitbleaching:stageadult                0.0633  0.3210   0.1973  0.8446 
traitsymbiont community:stageadult       0.2378  0.3573   0.6657  0.5094 
                                         ci.lb    ci.ub 
intrcpt                                -0.9204  -0.1926   ** 
traitnutrient content                  -0.0569   0.9071    . 
traitbleaching                         -0.4940   0.6557      
traitsymbiont community                -0.3908   0.7722      
traitsurvival                           0.1575   0.6748   ** 
stagejuvenile                          -1.1158  -0.1172    * 
stageadult                             -0.3209   0.4317      
traitbleaching:stagejuvenile           -0.8810   0.4733      
traitsymbiont community:stagejuvenile  -1.2357   1.6125      
traitsurvival:stagejuvenile             0.5245   1.1771  *** 
traitnutrient content:stageadult       -1.1812  -0.2446   ** 
traitbleaching:stageadult              -0.5854   0.7121      
traitsymbiont community:stageadult     -0.4843   0.9599      

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Plot the coefficients with different model datasets:

```r
dat_shg.out1 <- update(final.model.C, .~., data=dat_shg[-cdist.i[1],])
dat_shg.out2 <- update(final.model.C, .~., data=dat_shg[-cdist.i[2],])
dat_shg.out3 <- update(final.model.C, .~., data=dat_shg[-cdist.i[3],])
dat_shg.out4 <- update(final.model.C, .~., data=dat_shg[-cdist.i[4],])
dat_shg.out5 <- update(final.model.C, .~., data=dat_shg[-cdist.i[5],])
dat_shg.out.all <- update(final.model.C, .~., data=dat_shg[-cdist.i[c(1,2,3,5)],])


rbind(
	data.frame(
		Model = rep(" Original model fit"),
		Parameter = names(coef(final.model.C)), 
		Estimate = coef(final.model.C),
		SE = final.model.C$se),
	data.frame(
		Model = rep(paste0(" without Kenkel et al. 2015 - juvenile survival\nCook's dist = ",round(cdist[cdist.i[1]],1))),
		Parameter = names(coef(dat_shg.out1)), 
		Estimate = coef(dat_shg.out1),
		SE = dat_shg.out1$se),
	data.frame(
		Model = rep(paste0(" without Quigley et al. 2020 - juvenile survival\nCook's dist = ",round(cdist[cdist.i[2]],1))),
		Parameter = names(coef(dat_shg.out2)), 
		Estimate = coef(dat_shg.out2),
		SE = dat_shg.out2$se),
	data.frame(
		Model = rep(paste0(" without Manzello et al. 2019 - adult bleaching\nCook's dist = ",round(cdist[cdist.i[3]],1))),
		Parameter = names(coef(dat_shg.out3)), 
		Estimate = coef(dat_shg.out3),
		SE = dat_shg.out3$se),
	data.frame(
		Model = rep(paste0(" without Quigley et al. 2020 - juvenile bleaching\nCook's dist = ",round(cdist[cdist.i[4]],1))),
		Parameter = names(coef(dat_shg.out4)), 
		Estimate = coef(dat_shg.out4),
		SE = dat_shg.out4$se),
	data.frame(
		Model = rep(paste0(" without Császár et al. 2010 - adult growth\nCook's dist = ",round(cdist[cdist.i[5]],1))),
		Parameter = names(coef(dat_shg.out5)), 
		Estimate = coef(dat_shg.out5),
		SE = dat_shg.out5$se),
	data.frame(
		Model = rep(paste0("without estimates with Cook's distance > 5")),
		Parameter = names(coef(dat_shg.out.all)), 
		Estimate = coef(dat_shg.out.all),
		SE = dat_shg.out.all$se)) %>%
	mutate(Parameter = factor(Parameter),
		   Parameter = fct_relevel(Parameter, "intrcpt")) %>%
	ggplot(aes(x=Estimate, y=Parameter, color=Model)) +
	geom_vline(xintercept=0, linetype="dashed") +
	geom_pointrange(aes(xmin = Estimate - SE, xmax = Estimate + SE), position = position_dodge(width=0.3)) + 
	labs(x="Parameter estimate (log[X+0.2] scale) ± SE") +
	theme(legend.position="top", legend.direction="vertical")
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-39-1.png)<!-- -->


#### Residuals and simulated data


```r
plot(fitted(final.model.C), rstandard(final.model.C)$z)
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-40-1.png)<!-- -->

```r
plot(residuals(final.model.C, type="pearson"))
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-40-2.png)<!-- -->

```r
# Not amazing, not terrible

# Simulate new data given the parameter values of the data and same variance structure, plot on real data to see if it is similar
xsim <- simulate(final.model.C, 100) %>%
	mutate(true.val = dat_shg$val.log,
		   weighting = 1/dat_shg$sv.log) %>%
	pivot_longer(cols=-c(true.val, weighting)) %>%
	select(name, true.val, sim.val = value, weighting)
xsim %>%
	ggplot(aes(x=true.val, y=sim.val)) +
	geom_point(aes(col=weighting), alpha=0.05) +
	geom_abline(slope=1, intercept=0, linetype="dashed") + 
	scale_colour_viridis_c() +
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	scale_y_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-40-3.png)<!-- -->

```r
# reruns model, based around the same variance/covariance structures but with simulated data.
xsim %>% 
	ggplot(aes(y=weighting)) + 
	geom_point(aes(x=sim.val), col="red", size=2, alpha=0.05) +
	geom_point(data= filter(xsim, name == "sim_1"),
			   aes(x=true.val),
			   col="blue", size=2, alpha=0.9) +	
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	labs(x="Heritability estimate", y=expression("Weighting (1/SE"^2*")"))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-40-4.png)<!-- -->

Seems that the model is in general simulating the true data ok. Main problem is that values are simulating beyond the bounds of possible heritability.


<!-- ```{r, eval=F, include=F} -->
<!-- # How to do post-hoc tests -->
<!-- final.model.C <- rma(yi=val.log ~ trait*stage + growth.form,  -->
<!-- 					 vi=sv.log, data=dat_shg, method="REML", test="knha") -->
<!-- summary(final.model.C) -->

<!-- # Get contrasts:  -->
<!-- cont.mat <- rbind(c(0,0,0,0,0,0,0,0, 1,0, 0,0,0,0,0), # corymb vs. branch -->
<!-- 				  c(0,0,0,0,0,0,0,0, 0,1, 0,0,0,0,0), # corymb vs. mass -->
<!-- 				  c(0,0,0,0,0,0,0,0, 1,-1,0,0,0,0,0)) # branch vs. mass -->
<!-- # anova(final.model.D, L=cont.mat) # unadjusted -->
<!-- summary(multcomp::glht(final.model.C, -->
<!-- 					   linfct=cont.mat,  -->
<!-- 					   test=Ftest())) -->
<!-- 					   # test = chisq())) -->
<!-- 					   # test=adjusted("none"))) -->
<!-- ``` -->

------------------------------------------------------------------------

## Part D-i: Sub-analysis of trait type x **temperature differential**

### Step 0: View 'complete' dataset


```r
dat_tm %>% pull(temp.diff) %>% hist() # some high values, use a sqrt-transform to make more normally distributed!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

```r
dat_tm %>% pull(temp.diff) %>% sqrt() %>% hist() # better!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-41-2.png)<!-- -->

```r
dat_tm %>% table.plot("stage", "trait", .) # no other interactions possible!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-41-3.png)<!-- -->

```r
dat_tm %>% table.plot("h2", "trait", .) # no other interactions possible!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-41-4.png)<!-- -->

```r
dat_tm %>% table.plot("growth.form", "trait", .) # no other interactions
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-41-5.png)<!-- -->

Life stage, heritability type, and growth form are all signular for the photochemistry and immune response categories (all entirely overlapping for these trait types)
### Step 1: Define beyond optimal model

Beyond optimal is simply trait x temp differential:
this sub-analysis. $$ logit(h^2) \sim N(\mu, \tau^2) \\
\mu = trait_k \times sqrt(temp\;diff) + stage + h^2type + growth\;form \,  + \, \epsilon_{i} $$

### Step 2: Select optimal random effects


```r
m.study.nested <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
						 V=sv.log, random = ~1|study/est.id, 
						 data=dat_tm, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.spp.nested <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
						 V=sv.log, random = ~1|species/est.id, 
						 data=dat_tm, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.only.est <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm, method="REML", tdist=T,
					  sigma2 = c(0, NA)) # same as random = ~1|est.id
m.only.study <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|study
m.only.spp <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
					  V=sv.log, random = ~1|species/est.id, 
					  data=dat_tm, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|species
m.no.reff <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm, method="REML", tdist=T,
					  sigma2 = c(0, 0)) # no random effects


TableS9 <- myAIC(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff) %>%
	mutate(model = fct_recode(model,
		"Estimate ID nested in study ID" = "m.study.nested",
		"Estimate ID nested in species" = "m.spp.nested",
		"Study ID only" = "m.only.study",
		"Species only" = "m.only.spp",
		"Estimate ID only" = "m.only.est",
		"No random effect" = "m.no.reff"))
write.csv(TableS9, file="Suppl Tables/TableS9.csv", row.names=F)
TableS9
```

<div class="kable-table">

|model                          | df|     AICc|      ΔAICc|
|:------------------------------|--:|--------:|----------:|
|Estimate ID nested in study ID | 18| 106.2510|   0.000000|
|Estimate ID only               | 18| 109.9066|   3.655575|
|Estimate ID nested in species  | 18| 115.7686|   9.517644|
|Study ID only                  | 18| 122.1959|  15.944893|
|No random effect               | 18| 212.4151| 106.164141|
|Species only                   | 18| 212.4499| 106.198954|

</div>

```r
rm(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff)
```

The nested estimate ID within study ID random effects structure is preferred.
The next closest is estimate ID-only, with ΔAICc = 3.69.

### Step 3: Select optimal fixed effects


```r
mod_names <- ls(pattern = "m.fixed.")
do.call("rm", as.list(mod_names))

m.fixed.txm.s.h.g <- 
	rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.s.h <- 
	rma.mv(yi=val.log ~ trait*temp.s + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.s.g <- 
	rma.mv(yi=val.log ~ trait*temp.s + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.h.g <- 
	rma.mv(yi=val.log ~ trait*temp.s + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.s <- 
	rma.mv(yi=val.log ~ trait*temp.s + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.h <- 
	rma.mv(yi=val.log ~ trait*temp.s + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.g <- 
	rma.mv(yi=val.log ~ trait*temp.s + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm <- 
	rma.mv(yi=val.log ~ trait*temp.s,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))



m.fixed.t.m.s.h.g <- 
	rma.mv(yi=val.log ~ trait + temp.s + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.s.h <- 
	rma.mv(yi=val.log ~ trait + temp.s + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.s.g <- 
	rma.mv(yi=val.log ~ trait + temp.s + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.h.g <- 
	rma.mv(yi=val.log ~ trait + temp.s + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.s.h.g <- 
	rma.mv(yi=val.log ~ trait + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.s.h.g <- 
	rma.mv(yi=val.log ~ temp.s + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))


m.fixed.t.m.s <- 
	rma.mv(yi=val.log ~ trait + temp.s + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.h <- 
	rma.mv(yi=val.log ~ trait + temp.s + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.g <- 
	rma.mv(yi=val.log ~ trait + temp.s + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.s.h <- 
	rma.mv(yi=val.log ~ trait + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.s.g <- 
	rma.mv(yi=val.log ~ trait + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.h.g <- 
	rma.mv(yi=val.log ~ trait + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.m.s.h <- 
	rma.mv(yi=val.log ~ temp.s + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.s.g <- 
	rma.mv(yi=val.log ~ temp.s + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.h.g <- 
	rma.mv(yi=val.log ~ temp.s + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.s.h.g <- 
	rma.mv(yi=val.log ~ stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.t.m <- 
	rma.mv(yi=val.log ~ trait + temp.s,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.s <- 
	rma.mv(yi=val.log ~ trait + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.h <- 
	rma.mv(yi=val.log ~ trait + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.g <- 
	rma.mv(yi=val.log ~ trait + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.m.s <- 
	rma.mv(yi=val.log ~ temp.s + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.h <- 
	rma.mv(yi=val.log ~ temp.s + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.g <- 
	rma.mv(yi=val.log ~ temp.s + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.s.h <- 
	rma.mv(yi=val.log ~ stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.s.g <- 
	rma.mv(yi=val.log ~ stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.h.g <- 
	rma.mv(yi=val.log ~ h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))


m.fixed.t <- 
	rma.mv(yi=val.log ~ trait,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m <- 
	rma.mv(yi=val.log ~ temp.s,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.s <- 
	rma.mv(yi=val.log ~ stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.h <- 
	rma.mv(yi=val.log ~ h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.g <- 
	rma.mv(yi=val.log ~ growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.null <- 
	rma.mv(yi=val.log ~ 1,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

mod_names <- ls(pattern = "m.fixed.")
TableS10 <- head(do.call("myAIC", as.list(c(mod_names, getmodel=T)))) %>%
	mutate(model = fct_recode(model,
		"trait + heritability type" = "m.fixed.t.h",
		"trait" = "m.fixed.t",
		"trait x temp diff" = "m.fixed.txm",
		"trait + heritability type + temp diff" = "m.fixed.t.m.h",
		"trait + temp diff" = "m.fixed.t.m",
		"trait x temp diff + heritability type" = "m.fixed.txm.h"))
write.csv(TableS10, file="Suppl Tables/TableS10.csv", row.names=F)
TableS10
```

<div class="kable-table">

|model                                 | df|     AICc|     ΔAICc|
|:-------------------------------------|--:|--------:|---------:|
|trait + heritability type             |  6| 69.37282| 0.0000000|
|trait                                 |  5| 69.86893| 0.4961053|
|trait x temp diff                     | 11| 71.25972| 1.8869012|
|trait + heritability type + temp diff |  7| 71.68958| 2.3167581|
|trait + temp diff                     |  6| 72.36935| 2.9965266|
|trait x temp diff + heritability type | 12| 72.51076| 3.1379393|

</div>

### Step 4: Fit final model and alternative model including temperature


```r
final.model.Di <- rma.mv(yi=val.log ~ trait + h2, 
					  V=sv.log, random = ~1|study/est.id,
					  data=dat_tm, method="REML", tdist=T, sigma2 = c(NA, NA))
final.model.Di
```

```

Multivariate Meta-Analysis Model (k = 70; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0722  0.2688     12     no         study 
sigma^2.2  0.0616  0.2482     70     no  study/est.id 

Test for Residual Heterogeneity:
QE(df = 63) = 343.8923, p-val < .0001

Test of Moderators (coefficients 2:7):
F(df1 = 6, df2 = 63) = 8.0903, p-val < .0001

Model Results:

                       estimate      se     tval    pval    ci.lb    ci.ub 
intrcpt                 -0.7526  0.1600  -4.7042  <.0001  -1.0723  -0.4329  *** 
traitgrowth              0.0484  0.1488   0.3253  0.7460  -0.2490   0.3458      
traitnutrient content    0.1898  0.1764   1.0760  0.2860  -0.1627   0.5422      
traitbleaching           0.1699  0.2035   0.8350  0.4069  -0.2367   0.5765      
traitimmune response     0.5845  0.1969   2.9681  0.0042   0.1910   0.9780   ** 
traitsurvival            0.8209  0.1564   5.2501  <.0001   0.5084   1.1334  *** 
h2narrow                -0.3904  0.2275  -1.7159  0.0911  -0.8451   0.0643    . 

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
save(final.model.Di, file="Models/final.model.Di.Rdata")

mlm.variance.distribution(final.model.Di)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-44-1.png)<!-- --><div class="kable-table">

|        | % of total variance|I2    |
|:-------|-------------------:|:-----|
|Level 1 |            10.62483|---   |
|Level 2 |            41.13839|41.14 |
|Level 3 |            48.23679|48.24 |

</div>

```r
final.model.Di$k - final.model.Di$p
```

```
[1] 63
```

```r
# Temp x trait type model:
tempxtrait.model <- rma.mv(yi=val.log ~ trait * temp.s, 
					  V=sv.log, random = ~1|study/est.id,
					  data=dat_tm, method="REML", tdist=T, sigma2 = c(NA, NA))
tempxtrait.model
```

```

Multivariate Meta-Analysis Model (k = 70; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0992  0.3149     12     no         study 
sigma^2.2  0.0405  0.2013     70     no  study/est.id 

Test for Residual Heterogeneity:
QE(df = 58) = 257.6998, p-val < .0001

Test of Moderators (coefficients 2:12):
F(df1 = 11, df2 = 58) = 7.2539, p-val < .0001

Model Results:

                              estimate      se     tval    pval    ci.lb 
intrcpt                        -0.8284  0.2222  -3.7288  0.0004  -1.2731 
traitgrowth                     0.0443  0.2268   0.1954  0.8457  -0.4097 
traitnutrient content           0.1502  0.2405   0.6245  0.5347  -0.3313 
traitbleaching                  0.4375  0.2995   1.4608  0.1495  -0.1620 
traitimmune response            0.1466  0.2649   0.5533  0.5822  -0.3837 
traitsurvival                   0.8591  0.2187   3.9288  0.0002   0.4214 
temp.s                         -0.0079  0.1345  -0.0585  0.9535  -0.2771 
traitgrowth:temp.s             -0.0151  0.1586  -0.0951  0.9246  -0.3326 
traitnutrient content:temp.s    0.0526  0.1717   0.3062  0.7606  -0.2911 
traitbleaching:temp.s          -0.2549  0.1976  -1.2902  0.2021  -0.6505 
traitimmune response:temp.s     0.5337  0.2023   2.6382  0.0107   0.1288 
traitsurvival:temp.s           -0.0296  0.1519  -0.1952  0.8459  -0.3336 
                                ci.ub 
intrcpt                       -0.3837  *** 
traitgrowth                    0.4983      
traitnutrient content          0.6317      
traitbleaching                 1.0371      
traitimmune response           0.6768      
traitsurvival                  1.2969  *** 
temp.s                         0.2613      
traitgrowth:temp.s             0.3024      
traitnutrient content:temp.s   0.3962      
traitbleaching:temp.s          0.1406      
traitimmune response:temp.s    0.9387    * 
traitsurvival:temp.s           0.2743      

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# save(tempxtrait.model, file="Models/final.model.Di.Rdata")

mlm.variance.distribution(tempxtrait.model)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-44-2.png)<!-- --><div class="kable-table">

|        | % of total variance|I2    |
|:-------|-------------------:|:-----|
|Level 1 |            10.22495|---   |
|Level 2 |            26.04400|26.04 |
|Level 3 |            63.73104|63.73 |

</div>

```r
tempxtrait.model$k - tempxtrait.model$p
```

```
[1] 58
```


### Step 5: Check final model diagnostics (for alternative model)

#### Funnel plot


```r
x <- resid(tempxtrait.model)
hist(x) # looks approx. normal, except for one strong outlier at -1.5
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-45-1.png)<!-- -->

```r
par(mfrow=c(2,1))
funnel(tempxtrait.model)
update(tempxtrait.model, .~., data=dat_tm[x!=0,]) %>%
	funnel(., yaxis="vinv")
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-45-2.png)<!-- -->

```r
par(mfrow=c(1,1))
```

#### Fail-safe number


```r
fsn(yi=val.log, vi=sv.log, data=dat_tm, type="Rosenberg") # N=1669 to get no significance (but assumes fixed effects)
```

```

Fail-safe N Calculation Using the Rosenberg Approach

Average Effect Size:        -0.1361
Observed Significance Level: <.0001
Target Significance Level:   0.05

Fail-safe N: 1669
```

```r
# Number of studies x5 + 10:
length(unique(dat_tm$study))*5 + 10 # 1669 >> 70
```

```
[1] 70
```


#### Cook's distance


```r
cdist <- cooks.distance(tempxtrait.model, reestimate=FALSE)
plot(cdist)
cdist.i <- which(cdist >2)
points(cdist.i, cdist[cdist.i], col="red", pch=16)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-47-1.png)<!-- -->

```r
data[cdist.i,] %>% select(study, species, val, trait, stage, h2) %>%
	cbind(., cooks = round(cdist[cdist.i], 2))
```

<div class="kable-table">

|   |study               |species             |  val|trait           |stage    |h2     | cooks|
|:--|:-------------------|:-------------------|----:|:---------------|:--------|:------|-----:|
|3  |Quigley et al. 2020 |Acropora spathulata | 0.93|survival        |juvenile |narrow |  2.67|
|16 |Wright et al. 2019  |Acropora millepora  | 0.92|immune response |adult    |broad  |  7.64|

</div>

```r
update(tempxtrait.model, .~., data=dat_tm[-cdist.i,]) # effect of juv:growth n.s.
```

```

Multivariate Meta-Analysis Model (k = 68; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.1472  0.3836     12     no         study 
sigma^2.2  0.0089  0.0942     68     no  study/est.id 

Test for Residual Heterogeneity:
QE(df = 56) = 247.6234, p-val < .0001

Test of Moderators (coefficients 2:12):
F(df1 = 11, df2 = 56) = 7.9160, p-val < .0001

Model Results:

                              estimate      se     tval    pval    ci.lb 
intrcpt                        -0.8029  0.1882  -4.2673  <.0001  -1.1799 
traitgrowth                     0.0827  0.1717   0.4816  0.6319  -0.2613 
traitnutrient content           0.1550  0.1755   0.8834  0.3808  -0.1965 
traitbleaching                  0.4262  0.2343   1.8190  0.0743  -0.0432 
traitimmune response            0.1522  0.1972   0.7721  0.4433  -0.2428 
traitsurvival                   0.8945  0.1580   5.6622  <.0001   0.5780 
temp.s                         -0.0349  0.0989  -0.3532  0.7252  -0.2332 
traitgrowth:temp.s              0.0204  0.1188   0.1714  0.8645  -0.2176 
traitnutrient content:temp.s    0.0607  0.1285   0.4724  0.6385  -0.1968 
traitbleaching:temp.s          -0.1523  0.1586  -0.9606  0.3409  -0.4700 
traitimmune response:temp.s     0.2668  0.2114   1.2622  0.2121  -0.1567 
traitsurvival:temp.s           -0.0822  0.1171  -0.7022  0.4855  -0.3168 
                                ci.ub 
intrcpt                       -0.4260  *** 
traitgrowth                    0.4267      
traitnutrient content          0.5066      
traitbleaching                 0.8956    . 
traitimmune response           0.5473      
traitsurvival                  1.2109  *** 
temp.s                         0.1633      
traitgrowth:temp.s             0.2583      
traitnutrient content:temp.s   0.3182      
traitbleaching:temp.s          0.1653      
traitimmune response:temp.s    0.6903      
traitsurvival:temp.s           0.1523      

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
update(tempxtrait.model, .~., data=dat_tm[-cdist.i[1],]) # no real difference
```

```

Multivariate Meta-Analysis Model (k = 69; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.1437  0.3790     12     no         study 
sigma^2.2  0.0134  0.1156     69     no  study/est.id 

Test for Residual Heterogeneity:
QE(df = 57) = 254.1203, p-val < .0001

Test of Moderators (coefficients 2:12):
F(df1 = 11, df2 = 57) = 11.2090, p-val < .0001

Model Results:

                              estimate      se     tval    pval    ci.lb 
intrcpt                        -0.8102  0.1944  -4.1675  0.0001  -1.1996 
traitgrowth                     0.0859  0.1812   0.4740  0.6373  -0.2769 
traitnutrient content           0.1535  0.1862   0.8247  0.4130  -0.2192 
traitbleaching                  0.4308  0.2448   1.7603  0.0837  -0.0593 
traitimmune response            0.1509  0.2081   0.7250  0.4714  -0.2659 
traitsurvival                   0.8820  0.1683   5.2414  <.0001   0.5450 
temp.s                         -0.0273  0.1052  -0.2599  0.7958  -0.2379 
traitgrowth:temp.s              0.0148  0.1257   0.1180  0.9065  -0.2369 
traitnutrient content:temp.s    0.0568  0.1359   0.4179  0.6776  -0.2153 
traitbleaching:temp.s          -0.1600  0.1651  -0.9691  0.3366  -0.4906 
traitimmune response:temp.s     0.6038  0.1527   3.9550  0.0002   0.2981 
traitsurvival:temp.s           -0.0827  0.1234  -0.6702  0.5054  -0.3299 
                                ci.ub 
intrcpt                       -0.4209  *** 
traitgrowth                    0.4487      
traitnutrient content          0.5263      
traitbleaching                 0.9210    . 
traitimmune response           0.5677      
traitsurvival                  1.2190  *** 
temp.s                         0.1832      
traitgrowth:temp.s             0.2666      
traitnutrient content:temp.s   0.3289      
traitbleaching:temp.s          0.1706      
traitimmune response:temp.s    0.9096  *** 
traitsurvival:temp.s           0.1644      

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
update(tempxtrait.model, .~., data=dat_tm[-cdist.i[2],]) # no real difference
```

```

Multivariate Meta-Analysis Model (k = 69; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0994  0.3153     12     no         study 
sigma^2.2  0.0371  0.1926     69     no  study/est.id 

Test for Residual Heterogeneity:
QE(df = 57) = 251.2029, p-val < .0001

Test of Moderators (coefficients 2:12):
F(df1 = 11, df2 = 57) = 6.1569, p-val < .0001

Model Results:

                              estimate      se     tval    pval    ci.lb 
intrcpt                        -0.8252  0.2178  -3.7895  0.0004  -1.2612 
traitgrowth                     0.0415  0.2215   0.1872  0.8521  -0.4021 
traitnutrient content           0.1507  0.2344   0.6432  0.5227  -0.3186 
traitbleaching                  0.4355  0.2932   1.4854  0.1430  -0.1516 
traitimmune response            0.1471  0.2584   0.5693  0.5714  -0.3704 
traitsurvival                   0.8631  0.2130   4.0518  0.0002   0.4365 
temp.s                         -0.0115  0.1312  -0.0877  0.9304  -0.2742 
traitgrowth:temp.s             -0.0113  0.1549  -0.0730  0.9421  -0.3216 
traitnutrient content:temp.s    0.0548  0.1676   0.3271  0.7448  -0.2809 
traitbleaching:temp.s          -0.2545  0.1937  -1.3137  0.1942  -0.6424 
traitimmune response:temp.s     0.2434  0.2573   0.9459  0.3482  -0.2719 
traitsurvival:temp.s           -0.0270  0.1484  -0.1822  0.8560  -0.3242 
                                ci.ub 
intrcpt                       -0.3891  *** 
traitgrowth                    0.4850      
traitnutrient content          0.6201      
traitbleaching                 1.0226      
traitimmune response           0.6646      
traitsurvival                  1.2896  *** 
temp.s                         0.2512      
traitgrowth:temp.s             0.2989      
traitnutrient content:temp.s   0.3905      
traitbleaching:temp.s          0.1334      
traitimmune response:temp.s    0.7587      
traitsurvival:temp.s           0.2701      

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

#### Residuals and simulated data


```r
plot(fitted(tempxtrait.model), rstandard(tempxtrait.model)$z)
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

```r
plot(residuals(tempxtrait.model, type="pearson"))
abline(h=0, lty=2)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-48-2.png)<!-- -->

```r
# Not amazing, not terrible

# Simulate new data given the parameter values of the data and same variance structure, plot on real data to see if it is similar
xsim <- simulate(tempxtrait.model, 100) %>%
	mutate(true.val = dat_tm$val.log,
		   weighting = 1/dat_tm$sv.log) %>%
	pivot_longer(cols=-c(true.val, weighting)) %>%
	select(name, true.val, sim.val = value, weighting)
xsim %>%
	ggplot(aes(x=true.val, y=sim.val)) +
	geom_point(aes(col=weighting), alpha=0.05) +
	geom_abline(slope=1, intercept=0, linetype="dashed") + 
	scale_colour_viridis_c() +
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	scale_y_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-48-3.png)<!-- -->

```r
# reruns model, based around the same variance/covariance structures but with simulated data.
xsim %>% 
	ggplot(aes(y=weighting)) + 
	geom_point(aes(x=sim.val), col="red", size=2, alpha=0.05) +
	geom_point(data= filter(xsim, name == "sim_1"),
			   aes(x=true.val),
			   col="blue", size=2, alpha=0.9) +	
	scale_x_continuous(
		breaks=log( c(0,0.10,0.25,0.5,0.75,1)+0.2),
		labels =      c(0,0.10,0.25,0.5,0.75,1)) +
	labs(x="Heritability estimate", y=expression("Weighting (1/SE"^2*")"))
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-48-4.png)<!-- -->

------------------------------------------------------------------------


## Part D-ii: Sub-analysis of trait type x **binary temperature manipulation**

### Step 0: View 'complete' dataset


```r
dat_tm %>% table.plot("temp.manip", "trait", .)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-49-1.png)<!-- -->

```r
dat_tm %>% table.plot("temp.manip", "stage", .) # life stage x temp.manip possible
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-49-2.png)<!-- -->

```r
dat_tm %>% table.plot("temp.manip", "h2", .) # h2 x temp.manip
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-49-3.png)<!-- -->

```r
dat_tm %>% table.plot("temp.manip", "growth.form", .) # not possible for encrusting, columnar, remove!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-49-4.png)<!-- -->

```r
dat_tm %>% table.plot("stage", "trait", .) # no other interactions possible!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-49-5.png)<!-- -->

```r
dat_tm %>% table.plot("h2", "trait", .) # no other interactions possible!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-49-6.png)<!-- -->

```r
dat_tm %>% table.plot("growth.form", "trait", .) # no other interactions
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-49-7.png)<!-- -->

```r
# Life stage, heritability type, and growth form are all signular for the photochemistry and immune response categories (all entirely overlapping for these trait types)
```


### Step 1: Define beyond optimal model

Beyond optimal is simply trait x temp manipulation or temp magnitude for
this sub-analysis. $$ logit(h^2) \sim N(\mu, \tau^2) \\
\mu = trait_k \times temp\;manip + stage + h^2type + growth\;form \,  + \, \epsilon_{i} $$

### Step 2: Select optimal random effects


```r
## Compare random effects (fixed effects included)

m.study.nested <- rma.mv(yi=val.log ~ trait*temp.manip + stage + h2 + growth.form,
						 V=sv.log, random = ~1|study/est.id, 
						 data=dat_tm, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.spp.nested <- rma.mv(yi=val.log ~ trait*temp.manip + stage + h2 + growth.form,
						 V=sv.log, random = ~1|species/est.id, 
						 data=dat_tm, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.only.est <- rma.mv(yi=val.log ~ trait*temp.manip + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm, method="REML", tdist=T,
					  sigma2 = c(0, NA)) # same as random = ~1|est.id
m.only.study <- rma.mv(yi=val.log ~ trait*temp.manip + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|study
m.only.spp <- rma.mv(yi=val.log ~ trait*temp.manip + stage + h2 + growth.form,
					  V=sv.log, random = ~1|species/est.id, 
					  data=dat_tm, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|species
m.no.reff <- rma.mv(yi=val.log ~ trait*temp.manip + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm, method="REML", tdist=T,
					  sigma2 = c(0, 0)) # no random effects


TableS11 <- myAIC(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff) %>%
	mutate(model = fct_recode(model,
		"Estimate ID nested in study ID" = "m.study.nested",
		"Estimate ID nested in species" = "m.spp.nested",
		"Study ID only" = "m.only.study",
		"Species only" = "m.only.spp",
		"Estimate ID only" = "m.only.est",
		"No random effect" = "m.no.reff"))
write.csv(TableS11, file="Suppl Tables/TableS11.csv", row.names=F)
TableS11
```

<div class="kable-table">

|model                          | df|     AICc|     ΔAICc|
|:------------------------------|--:|--------:|---------:|
|Estimate ID nested in study ID | 18| 106.3283|  0.000000|
|Estimate ID only               | 18| 107.9884|  1.660187|
|Estimate ID nested in species  | 18| 113.8505|  7.522256|
|Study ID only                  | 18| 120.0003| 13.672022|
|Species only                   | 18| 178.0667| 71.738476|
|No random effect               | 18| 181.4887| 75.160427|

</div>

```r
rm(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff)
```

The nested estimate ID within study ID random effects structure is preferred.
The next closest is estimate ID-only, with ΔAICc = 1.66.

### Step 3: Select optimal fixed effects


```r
mod_names <- ls(pattern = "m.fixed.")
do.call("rm", as.list(mod_names))

m.fixed.txm.s.h.g <- 
	rma.mv(yi=val.log ~ trait*temp.manip + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.s.h <- 
	rma.mv(yi=val.log ~ trait*temp.manip + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.s.g <- 
	rma.mv(yi=val.log ~ trait*temp.manip + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.h.g <- 
	rma.mv(yi=val.log ~ trait*temp.manip + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.s <- 
	rma.mv(yi=val.log ~ trait*temp.manip + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.h <- 
	rma.mv(yi=val.log ~ trait*temp.manip + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm.g <- 
	rma.mv(yi=val.log ~ trait*temp.manip + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.txm <- 
	rma.mv(yi=val.log ~ trait*temp.manip,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))



m.fixed.t.m.s.h.g <- 
	rma.mv(yi=val.log ~ trait + temp.manip + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.s.h <- 
	rma.mv(yi=val.log ~ trait + temp.manip + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.s.g <- 
	rma.mv(yi=val.log ~ trait + temp.manip + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.h.g <- 
	rma.mv(yi=val.log ~ trait + temp.manip + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.s.h.g <- 
	rma.mv(yi=val.log ~ trait + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.s.h.g <- 
	rma.mv(yi=val.log ~ temp.manip + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))


m.fixed.t.m.s <- 
	rma.mv(yi=val.log ~ trait + temp.manip + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.h <- 
	rma.mv(yi=val.log ~ trait + temp.manip + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.m.g <- 
	rma.mv(yi=val.log ~ trait + temp.manip + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.s.h <- 
	rma.mv(yi=val.log ~ trait + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.s.g <- 
	rma.mv(yi=val.log ~ trait + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.h.g <- 
	rma.mv(yi=val.log ~ trait + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.m.s.h <- 
	rma.mv(yi=val.log ~ temp.manip + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.s.g <- 
	rma.mv(yi=val.log ~ temp.manip + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.h.g <- 
	rma.mv(yi=val.log ~ temp.manip + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.s.h.g <- 
	rma.mv(yi=val.log ~ stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.t.m <- 
	rma.mv(yi=val.log ~ trait + temp.manip,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.s <- 
	rma.mv(yi=val.log ~ trait + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.h <- 
	rma.mv(yi=val.log ~ trait + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.t.g <- 
	rma.mv(yi=val.log ~ trait + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.m.s <- 
	rma.mv(yi=val.log ~ temp.manip + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.h <- 
	rma.mv(yi=val.log ~ temp.manip + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m.g <- 
	rma.mv(yi=val.log ~ temp.manip + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.s.h <- 
	rma.mv(yi=val.log ~ stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.s.g <- 
	rma.mv(yi=val.log ~ stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.h.g <- 
	rma.mv(yi=val.log ~ h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))


m.fixed.t <- 
	rma.mv(yi=val.log ~ trait,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.m <- 
	rma.mv(yi=val.log ~ temp.manip,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.s <- 
	rma.mv(yi=val.log ~ stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.h <- 
	rma.mv(yi=val.log ~ h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))
m.fixed.g <- 
	rma.mv(yi=val.log ~ growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

m.fixed.null <- 
	rma.mv(yi=val.log ~ 1,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm, method="ML",
		   tdist=T, sigma = c(NA, NA))

mod_names <- ls(pattern = "m.fixed.")
TableS12 <- head(do.call("myAIC", as.list(c(mod_names, getmodel=T)))) %>%
	mutate(model = fct_recode(model,
		"trait + heritability type" = "m.fixed.t.h",
		"trait" = "m.fixed.t",
		"trait + heritability type + temp manip" = "m.fixed.t.m.h",
		"trait + temp manip" = "m.fixed.t.m",
		"trait x temp manip" = "m.fixed.txm",
		"trait + life stage" = "m.fixed.t.s"))
write.csv(TableS12, file="Suppl Tables/TableS12.csv", row.names=F)
TableS12
```

<div class="kable-table">

|model                                  | df|     AICc|     ΔAICc|
|:--------------------------------------|--:|--------:|---------:|
|trait + heritability type              |  6| 69.37282| 0.0000000|
|trait                                  |  5| 69.86893| 0.4961053|
|trait + heritability type + temp manip |  7| 70.17968| 0.8068565|
|trait + temp manip                     |  6| 71.22687| 1.8540525|
|trait x temp manip                     | 11| 71.67283| 2.3000087|
|trait + life stage                     |  7| 72.57154| 3.1987165|

</div>

A model of trait type + heritability type is preferred, though a model of trait (ΔAIC=0.5) and trait x temperature difference (ΔAIC=1.9) are also close.

### Step 4: Fit final model


```r
final.model.Dii <- rma.mv(yi=val.log ~ trait + h2, 
					  V=sv.log, random = ~1|study/est.id,
					  data=dat_tm, method="REML", tdist=T, sigma2 = c(NA, NA))
final.model.Dii
```

```

Multivariate Meta-Analysis Model (k = 70; method: REML)

Variance Components:

            estim    sqrt  nlvls  fixed        factor 
sigma^2.1  0.0722  0.2688     12     no         study 
sigma^2.2  0.0616  0.2482     70     no  study/est.id 

Test for Residual Heterogeneity:
QE(df = 63) = 343.8923, p-val < .0001

Test of Moderators (coefficients 2:7):
F(df1 = 6, df2 = 63) = 8.0903, p-val < .0001

Model Results:

                       estimate      se     tval    pval    ci.lb    ci.ub 
intrcpt                 -0.7526  0.1600  -4.7042  <.0001  -1.0723  -0.4329  *** 
traitgrowth              0.0484  0.1488   0.3253  0.7460  -0.2490   0.3458      
traitnutrient content    0.1898  0.1764   1.0760  0.2860  -0.1627   0.5422      
traitbleaching           0.1699  0.2035   0.8350  0.4069  -0.2367   0.5765      
traitimmune response     0.5845  0.1969   2.9681  0.0042   0.1910   0.9780   ** 
traitsurvival            0.8209  0.1564   5.2501  <.0001   0.5084   1.1334  *** 
h2narrow                -0.3904  0.2275  -1.7159  0.0911  -0.8451   0.0643    . 

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
save(final.model.Dii, file="Models/final.model.Dii.Rdata")

mlm.variance.distribution(final.model.Dii)
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-52-1.png)<!-- --><div class="kable-table">

|        | % of total variance|I2    |
|:-------|-------------------:|:-----|
|Level 1 |            10.62483|---   |
|Level 2 |            41.13839|41.14 |
|Level 3 |            48.23679|48.24 |

</div>

```r
final.model.Dii$k - final.model.Dii$p
```

```
[1] 63
```



------------------------------------------------------------------------

## Part D-iii: Sub-analysis of trait type x **temperature differential, with only non-zero values**

### Step 0: View 'complete' dataset


```r
dat_tm_nz <- dat_tm %>% filter(temp.diff >0) %>%
	mutate(trait = factor(trait), 
		   stage = factor(stage), 
		   h2 = factor(h2),
		   growth.form = factor(growth.form))
dat_tm_nz %>% pull(temp.diff) %>% hist() # some high values, use a sqrt-transform to make more normally distributed!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-53-1.png)<!-- -->

```r
dat_tm_nz %>% pull(temp.diff) %>% sqrt() %>% hist() # better!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-53-2.png)<!-- -->

```r
dat_tm_nz %>% table.plot("stage", "trait", .) # no other interactions possible!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-53-3.png)<!-- -->

```r
dat_tm_nz %>% table.plot("h2", "trait", .) # no other interactions possible!
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-53-4.png)<!-- -->

```r
dat_tm_nz %>% table.plot("growth.form", "trait", .) # no other interactions
```

![](/Users/Kribin/Documents/PhD Thesis/Heritability meta-analysis/Supplementary_Code_Documentation_B_Model_Selection_and_Diagnostics_2021-10-08_files/figure-html/unnamed-chunk-53-5.png)<!-- -->

```r
# Life stage, heritability type, and growth form are all singular for the photochemistry and immune response categories (all entirely overlapping for these trait types)
```


### Step 1: Define beyond optimal model

Beyond optimal is simply trait x temp differential:
this sub-analysis. $$ logit(h^2) \sim N(\mu, \tau^2) \\
\mu = trait_k \times sqrt(temp\;diff) + stage + h^2type + growth\;form \,  + \, \epsilon_{i} $$

### Step 2: Select optimal random effects


```r
m.study.nested <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
						 V=sv.log, random = ~1|study/est.id, 
						 data=dat_tm_nz, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.spp.nested <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
						 V=sv.log, random = ~1|species/est.id, 
						 data=dat_tm_nz, method="REML", tdist=T,
						 sigma2 = c(NA, NA))
m.only.est <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm_nz, method="REML", tdist=T,
					  sigma2 = c(0, NA)) # same as random = ~1|est.id
m.only.study <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm_nz, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|study
m.only.spp <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
					  V=sv.log, random = ~1|species/est.id, 
					  data=dat_tm_nz, method="REML", tdist=T,
					  sigma2 = c(NA, 0)) # same as random = ~1|species
m.no.reff <- rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
					  V=sv.log, random = ~1|study/est.id, 
					  data=dat_tm_nz, method="REML", tdist=T,
					  sigma2 = c(0, 0)) # no random effects


TableS13 <- myAIC(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff) %>%
	mutate(model = fct_recode(model,
		"Estimate ID nested in study ID" = "m.study.nested",
		"Estimate ID nested in species" = "m.spp.nested",
		"Study ID only" = "m.only.study",
		"Species only" = "m.only.spp",
		"Estimate ID only" = "m.only.est",
		"No random effect" = "m.no.reff"))
write.csv(TableS13, file="Suppl Tables/TableS13.csv", row.names=F)
TableS13
```

<div class="kable-table">

|model                          | df|     AICc|    ΔAICc|
|:------------------------------|--:|--------:|--------:|
|Estimate ID only               | 15| 120.4662|  0.00000|
|Study ID only                  | 15| 131.1002| 10.63401|
|Estimate ID nested in study ID | 15| 134.9756| 14.50947|
|Estimate ID nested in species  | 15| 137.2662| 16.80000|
|No random effect               | 15| 146.0231| 25.55699|
|Species only                   | 15| 150.1833| 29.71717|

</div>

```r
rm(m.study.nested, m.spp.nested, m.only.est, m.only.study, m.only.spp, m.no.reff)
```

The estimate ID-only random effects structure is preferred.
The next closest is estimate ID-only, with ΔAICc = 10.6.

### Step 3: Select optimal fixed effects


```r
mod_names <- ls(pattern = "m.fixed.")
do.call("rm", as.list(mod_names))

m.fixed.txm.s.h.g <- 
	rma.mv(yi=val.log ~ trait*temp.s + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.txm.s.h <- 
	rma.mv(yi=val.log ~ trait*temp.s + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.txm.s.g <- 
	rma.mv(yi=val.log ~ trait*temp.s + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.txm.h.g <- 
	rma.mv(yi=val.log ~ trait*temp.s + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.txm.s <- 
	rma.mv(yi=val.log ~ trait*temp.s + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.txm.h <- 
	rma.mv(yi=val.log ~ trait*temp.s + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.txm.g <- 
	rma.mv(yi=val.log ~ trait*temp.s + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.txm <- 
	rma.mv(yi=val.log ~ trait*temp.s,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))



m.fixed.t.m.s.h.g <- 
	rma.mv(yi=val.log ~ trait + temp.s + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.m.s.h <- 
	rma.mv(yi=val.log ~ trait + temp.s + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.m.s.g <- 
	rma.mv(yi=val.log ~ trait + temp.s + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.m.h.g <- 
	rma.mv(yi=val.log ~ trait + temp.s + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.s.h.g <- 
	rma.mv(yi=val.log ~ trait + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.m.s.h.g <- 
	rma.mv(yi=val.log ~ temp.s + stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))


m.fixed.t.m.s <- 
	rma.mv(yi=val.log ~ trait + temp.s + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.m.h <- 
	rma.mv(yi=val.log ~ trait + temp.s + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.m.g <- 
	rma.mv(yi=val.log ~ trait + temp.s + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.s.h <- 
	rma.mv(yi=val.log ~ trait + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.s.g <- 
	rma.mv(yi=val.log ~ trait + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.h.g <- 
	rma.mv(yi=val.log ~ trait + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))

m.fixed.m.s.h <- 
	rma.mv(yi=val.log ~ temp.s + stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.m.s.g <- 
	rma.mv(yi=val.log ~ temp.s + stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.m.h.g <- 
	rma.mv(yi=val.log ~ temp.s + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.s.h.g <- 
	rma.mv(yi=val.log ~ stage + h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))

m.fixed.t.m <- 
	rma.mv(yi=val.log ~ trait + temp.s,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.s <- 
	rma.mv(yi=val.log ~ trait + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.h <- 
	rma.mv(yi=val.log ~ trait + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.t.g <- 
	rma.mv(yi=val.log ~ trait + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))

m.fixed.m.s <- 
	rma.mv(yi=val.log ~ temp.s + stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.m.h <- 
	rma.mv(yi=val.log ~ temp.s + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.m.g <- 
	rma.mv(yi=val.log ~ temp.s + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))

m.fixed.s.h <- 
	rma.mv(yi=val.log ~ stage + h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.s.g <- 
	rma.mv(yi=val.log ~ stage + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.h.g <- 
	rma.mv(yi=val.log ~ h2 + growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))


m.fixed.t <- 
	rma.mv(yi=val.log ~ trait,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.m <- 
	rma.mv(yi=val.log ~ temp.s,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.s <- 
	rma.mv(yi=val.log ~ stage,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.h <- 
	rma.mv(yi=val.log ~ h2,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))
m.fixed.g <- 
	rma.mv(yi=val.log ~ growth.form,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))

m.fixed.null <- 
	rma.mv(yi=val.log ~ 1,
		   V=sv.log, random = ~1|study/est.id, data=dat_tm_nz, method="ML",
		   tdist=T, sigma =c(0, NA))

mod_names <- ls(pattern = "m.fixed.")
TableS14 <- head(do.call("myAIC", as.list(c(mod_names, getmodel=T)))) %>%
	mutate(model = fct_recode(model,
		"trait + life stage" = "m.fixed.t.s",
		"life stage" = "m.fixed.s",
		"trait + temp diff + life stage" = "m.fixed.t.m.s",
		"life stage + heritability type" = "m.fixed.s.h",
		"life stage + temp diff" = "m.fixed.m.s",
		"trait + life stage + heritability type" = "m.fixed.t.s.h"))
write.csv(TableS14, file="Suppl Tables/TableS14.csv", row.names=F)
TableS14
```

<div class="kable-table">

|model                                  | df|     AICc|     ΔAICc|
|:--------------------------------------|--:|--------:|---------:|
|trait + life stage                     |  7| 56.82920| 0.0000000|
|life stage                             |  2| 57.13176| 0.3025557|
|trait + temp diff + life stage         |  8| 57.39062| 0.5614156|
|life stage + heritability type         |  3| 59.35959| 2.5303886|
|life stage + temp diff                 |  3| 59.66963| 2.8404316|
|trait + life stage + heritability type |  8| 60.19922| 3.3700163|

</div>

### Step 4: Fit final model


```r
final.model.Diii <- rma(yi=val.log ~ trait + stage, vi=sv.log, 
					 data=dat_tm_nz, method="REML", test="knha")
final.model.Diii
```

```

Mixed-Effects Model (k = 44; tau^2 estimator: REML)

tau^2 (estimated amount of residual heterogeneity):     0.1017 (SE = 0.0350)
tau (square root of estimated tau^2 value):             0.3189
I^2 (residual heterogeneity / unaccounted variability): 75.32%
H^2 (unaccounted variability / sampling variability):   4.05
R^2 (amount of heterogeneity accounted for):            44.14%

Test for Residual Heterogeneity:
QE(df = 36) = 168.0695, p-val < .0001

Test of Moderators (coefficients 2:8):
F(df1 = 7, df2 = 36) = 4.6735, p-val = 0.0008

Model Results:

                       estimate      se     tval    pval    ci.lb    ci.ub 
intrcpt                 -0.3474  0.2144  -1.6206  0.1138  -0.7821   0.0874     
traitgrowth             -0.0745  0.2086  -0.3574  0.7229  -0.4975   0.3484     
traitnutrient content   -0.1846  0.2086  -0.8848  0.3821  -0.6077   0.2385     
traitbleaching          -0.2138  0.2411  -0.8869  0.3810  -0.7027   0.2751     
traitimmune response     0.4615  0.2687   1.7176  0.0945  -0.0834   1.0065   . 
traitsurvival            0.3768  0.1924   1.9587  0.0579  -0.0133   0.7670   . 
stagejuvenile           -0.6982  0.2096  -3.3317  0.0020  -1.1233  -0.2732  ** 
stageadult              -0.2565  0.1908  -1.3441  0.1873  -0.6435   0.1305     

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
save(final.model.Diii, file="Models/final.model.Diii.Rdata")

final.model.Diii$I
```

```
[1] 75.32359
```

```r
final.model.Diii$k - final.model.Diii$p
```

```
[1] 36
```



------------------------------------------------------------------------

R session info:


```r
sessionInfo()
```

```
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.4

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.0.1 metafor_2.4-0   Matrix_1.3-4    forcats_0.5.0  
 [5] stringr_1.4.0   dplyr_1.0.2     purrr_0.3.4     readr_1.3.1    
 [9] tidyr_1.1.2     tibble_3.0.3    ggplot2_3.3.5   tidyverse_1.3.0
[13] readxl_1.3.1   

loaded via a namespace (and not attached):
 [1] tidyselect_1.1.0  xfun_0.17         haven_2.3.1       lattice_0.20-41  
 [5] colorspace_1.4-1  vctrs_0.3.4       generics_0.0.2    viridisLite_0.3.0
 [9] htmltools_0.5.0   yaml_2.2.1        blob_1.2.1        rlang_0.4.11     
[13] pillar_1.4.6      glue_1.4.2        withr_2.2.0       DBI_1.1.0        
[17] dbplyr_1.4.4      modelr_0.1.8      lifecycle_0.2.0   munsell_0.5.0    
[21] gtable_0.3.0      cellranger_1.1.0  rvest_0.3.6       evaluate_0.14    
[25] labeling_0.3      knitr_1.29        fansi_0.4.1       highr_0.8        
[29] broom_0.7.0       Rcpp_1.0.5        scales_1.1.1      backports_1.1.10 
[33] jsonlite_1.7.2    farver_2.0.3      fs_1.5.0          hms_0.5.3        
[37] digest_0.6.28     stringi_1.5.3     grid_4.0.2        cli_2.0.2        
[41] tools_4.0.2       magrittr_2.0.1    crayon_1.3.4      pkgconfig_2.0.3  
[45] ellipsis_0.3.1    xml2_1.3.2        reprex_0.3.0      lubridate_1.7.9  
[49] assertthat_0.2.1  rmarkdown_2.5     httr_1.4.2        rstudioapi_0.11  
[53] R6_2.5.1          nlme_3.1-149      compiler_4.0.2   
```
