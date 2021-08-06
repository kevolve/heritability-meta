# Custom functions:

# load other custom functions:


gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}
# Mimics colours from the ggplot2 default colour wheel
# gg_color_hue(3)


# need a custom function to deal with range expansion bug for sqrt transformation:
mysqrt_trans <- function() {
	scales::trans_new("mysqrt", 
					  transform = base::sqrt,
					  inverse = function(x) ifelse(x<0, 0, x^2),
					  domain = c(0, Inf))
}


table.plot <- function(X = "trait", Y = "study", Data = data) {
	table.data <- Data[,c(X, Y)] %>% 
		table() %>% as_tibble %>%
		rename(X = all_of(X), Y = all_of(Y)) %>% 
		mutate(not.zero=case_when(n>0 ~ 1, TRUE ~ 0))
	max.val <- table.data %>% pull(n) %>% max()
	# col.vals <- if(any(table.data$not.zero == 0)) {
	# 	c("grey", scales::viridis_pal(option="C")(max.val-1))
	# 	} else {scales::viridis_pal(option="C")(max.val)}
	table.data %>% 
		mutate(n = ifelse(n==0, NA_real_, n)) %>%
		ggplot(aes(x = X, y = Y)) +
		geom_tile(aes(fill = n)) +
		scale_fill_viridis_c(na.value = "grey") +
		# scale_fill_continuous(values = col.vals) +
		# scale_fill_manual(values=color.vals, breaks=levels(ddf$zz)[seq(1, N, by=5)])
		# ?viridis::scale_fill_viridis(option="C") +
		geom_text(aes(label = replace_na(n, replace = 0), alpha=factor(not.zero, levels=c(1,0))), 
				  color = "black", fontface = "bold", size = 6) +
		scale_alpha_manual(values=c(1,0)) +
		guides(alpha=FALSE, fill=FALSE) +
		theme_bw() + 
		labs(x = "", y = "", 
			 title = paste0("n per ", X, " and ",Y,"\n","(N = ", sum(table.data$n), ")")) +
		theme(plot.title = element_text(hjust = 0.5, face="bold", 
										colour = "black", size = 16), 
			  plot.margin=unit(c(7,17,0,0), "pt"),
			  axis.text = element_text(face="bold", colour="black", size = 12),
			  axis.text.x=element_text(angle = 300, vjust=0.5, hjust=0))
}
# Plots matrix of cells based on completeness of data
## Add functionality: include number of independent traits/total for each study
## and similarly number of independent studies/total on opposite side of grid


summarise.coverage <- function(X = "trait", Y = "study", Data = data){
	Data <- Data %>% mutate(x = factor(eval(parse(text=X))), 
							y = factor(eval(parse(text=Y))))
	table_by_x <- Data %>% group_by(x, y) %>% count() %>% ungroup(x) %>%
		mutate(n=rep(1)) %>% count(wt=n, sort=TRUE) %>% as.data.frame
	colnames(table_by_x) <- c(Y, paste0("# ",X, " levels covered"))
	table_by_y <- Data %>% group_by(x, y) %>% count() %>% ungroup(y) %>%
		mutate(n=rep(1)) %>% count(wt=n, sort=TRUE) %>% as.data.frame
	colnames(table_by_y) <- c(X, paste0("# ",Y, " levels covered"))
	list(X = table_by_x, Y = table_by_y)
}
# Plots marginal counts/coverage between different factor levels
# summarise.coverage("trait", "study", data)


myAIC <- function(..., aic=F, aicc=T, waic=F, 
				  delta_aic=F, delta_aicc=T, pvals=F, 
				  getmodel=F, param.list = NULL){
	require(tidyverse)
	if(getmodel==T){
		mod.name <- as.character(list(...))
		mod.out <- list(...)
		mod.out <- lapply(mod.out, function(x) get(x))
	} else {
		mod.name <- eval(substitute(alist(...))) %>% as.character
		mod.out <- list(...)
	}
	
	# add additional optional arguments later...
	dfs <- mod.out %>% map(anova) %>% map_dbl("m")
	AICs <- mod.out %>% map_dbl(AIC)
	AICcs <- mod.out %>% map_dbl(AIC, correct=T)
	p <- mod.out %>% map(anova) %>% map_dbl("QMp")
	
	out <- data.frame(model=mod.name, df=dfs, AIC=AICs, AICc=AICcs, p=p) %>%
		mutate(ΔAIC = AIC - min(AIC),
				ΔAICc = AICc - min(AICc),
				p = format(p, scientific=F, digits=2),
				expAICc = exp(-0.5*ΔAICc),
				wAICc = expAICc/max(expAICc),
				expAICc = NULL) %>%
		relocate(ΔAIC, .after="AIC") %>%
		relocate(ΔAICc, .after="AICc") %>%
		rename("Overall p-val"=p) %>% 
		arrange(AICc)
	
	if(delta_aic==T) {out <- out %>% arrange(AIC)}
	
	if(!is.null(param.list)) {
		possible.terms <- paste("x ~",paste(param.list, collapse="*")) %>%
			map(as.formula) %>% map(terms.formula) %>%
			map(~attr(.x, "term.labels")) %>% unlist
		chars <- mod.out %>% map("formula.yi") %>% as.character()
		terms <- chars %>% map(as.formula) %>% map(terms.formula) %>%
			map(~attr(.x, "term.labels"))
		
		x <- matrix(nrow=0, ncol=length(possible.terms))
		for (i in 1:length(mod.name)){
			x <- rbind(x, possible.terms %in% terms[[i]] %>% as.numeric)
		}
		
		i.x <- colSums(x)>0
		x <- as.data.frame(x[,i.x]) * out$wAICc
		colnames(x) <- possible.terms[i.x]
		
		rel.importance <- x %>%	summarise_all(sum) %>%
			pivot_longer(everything()) %>%
			mutate(norm_val = value/max(value)) %>%
			arrange(-norm_val)
	}
	
	if(aic == F) {out <- out %>% select(-AIC)}
	if(aicc== F) {out <- out %>% select(-AICc)}
	if(delta_aic==F) {out <- out %>% select(-ΔAIC)}
	if(delta_aicc==F) {out <- out %>% select(-ΔAICc)}
	if(waic==F) {out <- out %>% select(-wAICc)}
	if(pvals==F) {out <- out %>% select(-`Overall p-val`)}
	
	if(!is.null(param.list)) {
		out <- list(AICtable=out, rel.importance=rel.importance)
	}
	out
}
# Usage:
# myAIC(m.h2, m.trait, m.growth.form, m.level, pvals=T)
# myAIC(m.h2, m.trait, m.growth.form, m.level, pvals=T,
# 	  param.list=c("trait", "h2", "growth.form", "stage", "level")) # not run
# myAIC("m.h2", "m.trait", "m.growth.form", "m.level", pvals=T, getmodel=T,
# 	  param.list=c("trait", "h2", "growth.form", "stage", "level"))
# do.call("myAIC", as.list(c(mod_names=mod_names, getmodel=T, param.list = list(param.list))))

make_modsel_table <- function(reff.tbl, feff.tbl, data.subset, beyond.opt.mod, n=1) {
	require(gt)
	require(tidyverse)
	bind_rows(reff.tbl, feff.tbl, .id = 'groupname') %>%
		mutate(groupname = fct_recode(groupname, "Top random effects"="1", "Top fixed effects"="2")) %>%
		group_by(groupname) %>%
		mutate(rank=1:n()) %>%
		gt(data=., groupname_col="groupname", rowname_col = "rank") %>%
		cols_label(model = "") %>%
		tab_header(title = md(paste0("**Table S",n,".** Model selection for random and fixed effect structures using ",data.subset,", ordered by lowest AICc."))) %>%
		tab_footnote(
			footnote = paste0("Models fit via REML using a 'beyond optimal' fixed effects structure of: ",beyond.opt.mod,"."),
			locations = cells_row_groups(groups="Top random effects")) %>%
		tab_footnote(
			footnote = paste("Models fit via ML using the top random effects structure, seen above."),
			locations = cells_row_groups(groups="Top fixed effects")) %>%
		fmt_number(columns = 4,	decimals = 1) %>%
		fmt_number(columns = 5,	decimals = 2) %>%
		tab_style(
			style = cell_text(weight="bold"),
			locations = cells_body(rows = (ΔAICc<2))) %>%
		tab_style(
			style = cell_text(weight="bold", style="italic"),
			locations = cells_body(columns = "model", rows = (ΔAICc==0))) %>%
		tab_style(
			style = cell_text(weight="bold"),
			locations = cells_stub(rows = c(1, nrow(reff.tbl)+1))) %>%
		tab_style(
			style = cell_text(align = "left", size="medium"),
			locations = cells_title(groups="title")) %>%
		tab_options(heading.border.bottom.color = "#989898",
					row_group.background.color = "#D3D3D3D3",
					row_group.border.top.color = "#D3D3D3D3",
					row_group.border.bottom.style = "none",
					stub.border.style = "none",
					footnotes.marks = "standard")
}
# Usage:
# x <- read.csv(file="GITIGNORE/Tables/TableS1A.csv") %>% head(3) 
# y <- read.csv(file="GITIGNORE/Tables/TableS1B.csv") %>% head(5)
# 
# make_modsel_table(x,y, n=1, data.subset = "the complete dataset",
# 		  beyond.opt.mod="trait + life stage + heritability type + growth form")


make_coef_table <- function(top.model, feff.tbl, data.subset, n=1, reff.footnote=NULL){
	
	require(gt)
	require(tidyverse)
	require(metafor)
	
	nb.models <- feff.tbl %>% 
		filter(ΔAICc < 2) %>% 
		select(model) %>% 
		mutate(model.f = 
			   	gsub(x=model, pattern="heritability type", replacement="h2") %>%
			   	gsub(x=., pattern="life stage", replacement="stage") %>%
			   	gsub(x=., pattern="growth form", replacement="growth.form") %>%
			   	gsub(x=., pattern="temp manip", replacement="temp.s") %>%
			   	gsub(x=., pattern="temp diff", replacement="temp.s") %>%
			   	gsub(x=., pattern="x", replacement="*") %>%
			   	paste0(". ~ ",.)) %>%
		mutate(rank = paste0("m",1:n()))
	nb.model.params <- nb.models %>%
		split(.$model.f) %>%
		map(~update(top.model, formula=.x$model.f)) %>%
		map_dfr(~coef(summary(.)), .id="model.f") %>%
		rownames_to_column("parameter") %>%
		mutate(parameter = 
			   	gsub(x=parameter, pattern="[^A-Za-z]", replacement="") %>%
			   	gsub(x=., pattern="trait|stage", replacement="")) %>%
		left_join(select(nb.models, model, model.f, rank)) %>%
		select(parameter, rank, model, estimate, se, ci.lb, ci.ub, tval, pval) %>%
		mutate(parameter = fct_relevel(parameter, unique(parameter))) %>%
		arrange(parameter, rank)
	
	coef.plot <- nb.model.params %>%
		ggplot(aes(x=estimate, y=parameter, color=model)) +
		geom_vline(xintercept=0, linetype="dashed") +
		geom_pointrange(aes(xmin = estimate - se, xmax = estimate + se), position = position_dodge(width=0.3)) + 
		labs(x="Parameter estimate (log[X+0.2] scale) ± SE", y=NULL) +
		theme(legend.position="top", legend.direction="vertical")
	
	if(is.null(reff.footnote)){
		reff.footnote <- "Models fit via REML using the optimal random effects structure selected previously."}
	
	out.table <- nb.model.params %>%
		mutate(model = fct_relevel(model, nb.models$model),
			   model = paste0("m",as.numeric(model))) %>%
		pivot_longer(estimate:pval) %>%
		unite("columns", name, model) %>% 
		pivot_wider(id_cols = parameter, names_from="columns", values_from = "value") %>%
		relocate(parameter, starts_with("estimate"), starts_with("se")) %>%
		relocate(starts_with("tval"), starts_with("pval"), .after=last_col()) %>%
		gt(rowname_col = "parameter") %>%
		tab_header(title = md(paste0("**Table S",n,".** Model coefficients (estimated on the log[h2+0.2]-scale) for the top ",
									 ifelse(nrow(nb.models)==1, "model (no other models within 2 ΔAICc)", 
									 	   paste(nrow(nb.models), "models within 2 ΔAICc")),
									 ", using ", data.subset,"."))) %>%
		fmt_number(columns=everything(), decimals = 2) %>%
		fmt_missing(columns=everything(), missing_text = "") %>%
		cols_merge(columns=vars(estimate_m1, se_m1, ci.lb_m1, ci.ub_m1), pattern="{1}±{2}<br>({3}, {4})") %>%
		tab_style(locations = cells_body(columns="estimate_m1", rows = is.na(estimate_m1)), style = cell_text(size=0)) %>%
		cols_merge(columns=vars(tval_m1, pval_m1), pattern="<i>t</i> = {1}<br><i>p</i> = {2}") %>%
		tab_style(locations = cells_body(columns="tval_m1", rows = is.na(tval_m1)), style = cell_text(size=0)) %>%
		tab_spanner(label=paste(nb.models$model[1],"model"), columns=vars(estimate_m1, tval_m1)) %>%
		cols_label(estimate_m1 = "est ± SE\n(95%CI)", tval_m1 = "test statistic") %>%
		tab_options(heading.border.bottom.color = "#989898", stub.border.style = "none",
					footnotes.border.bottom.color="#D3D3D3", footnotes.marks = "standard") %>% 
		tab_footnote(footnote = reff.footnote, locations = cells_title("title"))
	
	
	
	
	
	if(nrow(nb.models)>=2) {
		out.table <- out.table %>% 
			cols_merge(columns=vars(estimate_m2, se_m2, ci.lb_m2, ci.ub_m2), pattern="{1}±{2}<br>({3}, {4})") %>%
			tab_style(locations = cells_body(columns="estimate_m2", rows = is.na(estimate_m2)), style = cell_text(size=0)) %>%
			cols_merge(columns=vars(tval_m2, pval_m2), pattern="<i>t</i> = {1}<br><i>p</i> = {2}") %>%
			tab_style(locations = cells_body(columns="tval_m2", rows = is.na(tval_m2)), style = cell_text(size=0)) %>%
			tab_spanner(label=paste(nb.models$model[2],"model"), columns=vars(estimate_m2, tval_m2)) %>%
			cols_label(estimate_m2 = "est ± SE\n(95%CI)", tval_m2 = "test statistic")
		
	}
	if(nrow(nb.models)>=3) {
		out.table <- out.table %>% 
			cols_merge(columns=vars(estimate_m3, se_m3, ci.lb_m3, ci.ub_m3), pattern="{1}±{2}<br>({3}, {4})") %>%
			tab_style(locations = cells_body(columns="estimate_m3", rows = is.na(estimate_m3)), style = cell_text(size=0)) %>%
			cols_merge(columns=vars(tval_m3, pval_m3), pattern="<i>t</i> = {1}<br><i>p</i> = {2}") %>%
			tab_style(locations = cells_body(columns="tval_m3", rows = is.na(tval_m3)), style = cell_text(size=0)) %>%
			tab_spanner(label=paste(nb.models$model[3],"model"), columns=vars(estimate_m3, tval_m3)) %>%
			cols_label(estimate_m3 = "est ± SE\n(95%CI)", tval_m3 = "test statistic")
	}
	if(nrow(nb.models)>=4) {
		out.table <- out.table %>% 
			cols_merge(columns=vars(estimate_m4, se_m4, ci.lb_m4, ci.ub_m4), pattern="{1}±{2}<br>({3}, {4})") %>%
			tab_style(locations = cells_body(columns="estimate_m4", rows = is.na(estimate_m4)), style = cell_text(size=0)) %>%
			cols_merge(columns=vars(tval_m4, pval_m4), pattern="<i>t</i> = {1}<br><i>p</i> = {2}") %>%
			tab_style(locations = cells_body(columns="tval_m4", rows = is.na(tval_m4)), style = cell_text(size=0)) %>%
			tab_spanner(label=paste(nb.models$model[4],"model"), columns=vars(estimate_m4, tval_m4)) %>%
			cols_label(estimate_m4 = "est ± SE\n(95%CI)", tval_m4 = "test statistic")
	}
	print(coef.plot)
	out.table
}
# Usage:
# y <- read.csv(file="GITIGNORE/Tables/TableS6B.csv") %>% head(5)
# make_coef_table(final.model, y, data.subset = "the complete dataset", n=8)

