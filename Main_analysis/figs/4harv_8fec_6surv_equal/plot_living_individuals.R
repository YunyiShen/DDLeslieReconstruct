generate_plot_data <- function(w, group){
  CI_low <- apply(w, 1, quantile, 0.025)
  CI_high <- apply(w, 1, quantile, 0.975)
  
  data.frame(
    group = group, year = 1:17+1991,
    mean = rowMeans(w),
    CI_low = CI_low,
    CI_high = CI_high
  )
}


living_matrix <- lapply(1:2000, function(i, living_mcmc){
  matrix(living_mcmc[i,], nrow = 11)
}, Chicago_RES$mcmc.objs$living.mcmc)

all_inid <- sapply(living_matrix, function(w){
  t(colSums(w))
})

males <- sapply(living_matrix, function(w){
  t(colSums(w[10:11,])/colSums(w))
})

females <- sapply(living_matrix, function(w){
  t(colSums(w[2:8,])/colSums(w))
})

fawns <- sapply(living_matrix, function(w){
  t(colSums(w[c(1,9),])/colSums(w))
})

fawns <- generate_plot_data(fawns, "fawn")
females <- generate_plot_data(females, "female")
males <- generate_plot_data(males, "male")

plot_data <- rbind(fawns, females, males)
write.csv(plot_data, "post_cull_population_age_structure.csv", row.names = F)

library(ggplot2)
ggplot(plot_data, aes(x = year, y = mean, shape = group, lty = group))+
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),size = .5, width = 0.3)+
  xlab("Year") + 
  ylab("Reconstructed post-cull\n age-sex structure")+
  theme_classic() 
ggsave("postcull_ratio.pdf", width = 6, height = 3.5, scale = .9)


