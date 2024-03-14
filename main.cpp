#include <stddef.h>
#include <stdio.h>
#include <random>
#include <vector>
#include <math.h>

struct RandomValues {
    double u;
    double v;
    double s;
};

RandomValues genRand(std::mt19937 &eg, std::uniform_real_distribution<> &udist) {
    double s = 0.0;
    double u = 0.0;
    double v = 0.0;

    while (s == 0.0 || s >= 1.0) {
        u = udist(eg);
        v = udist(eg);

        s = u * u + v * v;
    }

    return { u, v, s };
}

std::vector<double> normallyDistributedValues(size_t npairs) {
    std::random_device rd;
    std::mt19937 eg(rd());
    std::uniform_real_distribution<> udist(-1.0, 1.0);

    std::vector<double> normVals;
    normVals.reserve(npairs * 2);

    for (size_t i = 0; i < npairs; i++) {
        RandomValues r = genRand(eg, udist);

        double z1 = r.u * sqrt(-2.0 * log(r.s) / r.s);
        double z2 = r.v * sqrt(-2.0 * log(r.s) / r.s);

        normVals.push_back(z1);
        normVals.push_back(z2);
    }

    return normVals;
}

size_t bucketingFunc(double x, double lo, double interval) {
    double lowerBound = lo - interval / 2.0;

    double y = floor((x - lowerBound) / interval);
    return static_cast<size_t>(y);
}

std::vector<double> histogram(std::vector<double> &values, double lo, double hi, double interval, size_t total) {
    size_t numBuckets = bucketingFunc(hi, lo, interval) + 1;
    std::vector<int> buckets(numBuckets, 0);

    double halfInterval = interval / 2.0;
    for (double val : values) {
        if (val <= lo - halfInterval || val >= hi + halfInterval) continue;

        size_t buck = bucketingFunc(val, lo, interval);
        buckets[buck]++;
    }

    std::vector<double> hist;
    for (int buck : buckets) {
        hist.push_back(static_cast<double>(buck) / total / interval);
    }

    return hist;
}

// y = (1/sqrt(2pi))e^(-x^2/2)

double fitNormal(double lo, double interval, std::vector<double> &yvals) {
    double totalX = 0.0;
    double totalY = 0.0;
    size_t n = yvals.size();

    for (size_t i = 0; i < n; i++) {
        double x = lo + i * interval;
        totalX += exp(-x * x / 2.0);
        totalY += yvals[i];
    }

    double meanX = totalX / n;
    double meanY = totalY / n;

    double totalNumerator = 0.0;
    double totalDenominator = 0.0;

    for (size_t i = 0; i < n; i++) {
        double x = lo + i * interval;

        double xdev = exp(-x * x / 2.0) - meanX;
        double ydev = yvals[i] - meanY;

        totalNumerator += xdev * ydev;
        totalDenominator += xdev * xdev;
    }

    return totalNumerator / totalDenominator;
}

int main() {

    size_t npairs = 5000000;
    auto normVals = normallyDistributedValues(npairs);
    
    double lowestBucket = 0.0;
    double highestBucket = 5.0;
    double interval = 0.1;
    auto approxNDist = histogram(normVals, lowestBucket, highestBucket, interval, npairs * 2);

    // ~ 1 / sqrt(2pi)
    double gradient = fitNormal(lowestBucket, interval, approxNDist);
    
    double pi = 1.0 / (2.0 * gradient * gradient);

    printf("%f\n", pi);
}