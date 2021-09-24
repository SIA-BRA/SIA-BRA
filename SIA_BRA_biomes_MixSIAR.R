## SIA-BRA_biomes
# Diniz-Reis et al 2021
# GEB
# MixSiar

install.packages("MixSIAR")
install.packages("ggplot2")

library(MixSIAR)
library(mcmc)
library(ggplot2)

# Consumers: 'fish_consumer', 'mammal_consumer','bird_consumer'

mix.filename <- system.file("extdata", "consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename, 
                     iso_names=c("d13C"), 
                     factors=c("biome"), 
                     fac_random=c(TRUE), 
                     fac_nested=c(FALSE), 
                     cont_effects=NULL)

# Source: C3 and C4 plants
source.filename <- system.file("extdata", "source_bra.csv", package = "MixSIAR")

source <- load_source_data(filename=source.filename, source_factors="biome", 
                           conc_dep=FALSE, data_type="means", mix)

# Discrimination: zero = d13Cd and d15Nd
discr.filename <- system.file("extdata", "discriminatrion_bra.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)

plot_prior(alpha.prior=1,source)


# model
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)


jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior = 1)

output_options <- list(summary_save = TRUE,                 
                       summary_name = "summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = FALSE,
                       return_obj = TRUE)
output_JAGS(jags.1, mix, source, output_options)

# Plot contributions - each consumer

grafico<-g.post$fac1[[4]]+ scale_x_continuous(breaks = seq(0,1,0.1), limits=c(0, 1))+
  scale_colour_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9")) +
  scale_fill_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9"))+ 
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
      axis.line.y = element_line(size = 0.5, colour = "black"),
      axis.line = element_line(size=1, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x=element_text(size=12, color="black"),
      axis.text.y=element_text(size=12, color="black"),
      axis.title.x=element_text(size=18, face="bold"),
      axis.title.y=element_text(size=16, face="bold"),
      legend.title = element_blank(),
      legend.text = element_text(size=14, face="bold"))
grafico
ggsave("mixsiar_class_consumer.jpeg", units="cm", width=14, height=12, dpi=300)
dev.off()