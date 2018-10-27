undig <- structure(list(Percent = c(0.058, 0.108, 0.047, 0.088, 0.086, 
0.15, 0.091, 0.16), Method = c("Conservative Estimate", "Upper Bound", 
"Conservative Estimate", "Upper Bound", "Conservative Estimate", 
"Upper Bound", "Conservative Estimate", "Upper Bound"), Race = structure(c(1L, 
1L, 5L, 5L, 3L, 2L, 4L, 4L), .Label = c("All MSM", "Black MSM", 
"Black MSM", "Hispanic MSM", "White MSM"), class = "factor")), .Names = c("Percent", 
"Method", "Race"), row.names = c("1", "2", "3", "4", "5", "6", 
"7", "8"), class = "data.frame")

ggplot() +
	geom_hline(data=undig,colour = '#ff9999',size = 0.75,alpha = 0.5052,yintercept = 1.0) +
	geom_point(aes(x = Percent,y = Race,shape = Method),data=undig,size = 4.0) +
	scale_x_continuous(name = 'Undiagnosed %',labels = percent_format(),limits = c(0,.20)) +
	theme_bw() +
	theme(legend.position = 'bottom') +
	scale_y_discrete(name = ' ')