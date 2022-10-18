## FAP experiments

This section stores the FAP values obtained during experiments.

- Periods tested are: {1, 4, 7, 10, 13, 16, 19}
- Transit duration tested are: {0.5, 1, 2, 3, 4, 5, 6} hours.
- Depths at which FAPs are calculated are:

	```
	depths <- c(
    	seq(from = 0.005, to = 0.03, length.out=30),
    	seq(from = 0.031, to = 0.1, length.out=10)
	)
	```
