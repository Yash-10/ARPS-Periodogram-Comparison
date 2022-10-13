## FAP experiments

This section stores the FAP values obtained during experiments.

- Folders of type `transit_duration_xhr` use L=500, R=300, ntransits=10, significanceMode='max', noiseType=1, periods from {1, 3, 5, 7, 9, 11} days, and use depths:

```
depths <- c(
    seq(from = 0.005, to = 0.03, length.out=30),
    seq(from = 0.031, to = 0.1, length.out=10)
)
```
