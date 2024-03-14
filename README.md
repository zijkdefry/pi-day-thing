**Pi Day**

Boring (and probably unoriginal) Pi Day thing I did for 2024-03-14

```
$ ./piday
3.142407
```

Summary:
1. By using the [Box-Muller Transform](https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform),
two uniformly sampled random values can be used to obtain two other values that are normally distributed.

2. From this, the normal distrubution function can be approximated by sampling the above many times.

3. Since the normal distribution has the form `y = e^(-x^2 / 2) / 2sqrt(PI)`,
the graph of `y` against `e^(-x^2 / 2)` has gradient `1 / 2sqrt(PI)`.

4. Fitting our approximation of the normal distribution allows us to calculate the approximate
value of `PI` from the gradient of the best fit line.

\>3,000,000 samples are needed to get a consistent 2-decimal-place accuracy of `PI`