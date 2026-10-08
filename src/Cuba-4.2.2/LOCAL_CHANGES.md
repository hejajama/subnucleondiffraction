# Local changes to Cuba 4.2.2

This copy of Cuba differs from the release. When updating Cuba, check
whether the release fixes these issues, and reapply the changes otherwise.

## suave/Sample.c: region without a usable sample set

Suave adds the weight of a sample set to a region only if the set has at
least `nmin` points. After a split, a half-region can end up without such a
set, for example where the integrand is (practically) zero, and then its
weight sum is 0. Suave computed `sigsq = 1/weightsum` and
`avg = sigsq*avgsum`, which gives infinity and NaN. The NaN spreads into the
total, and through the fluctuations into the number of samples of the next
split, which gives NaN results, failed allocations or a "buffer overflow
detected" abort.

The change sets the result and the variance of such a region to 0. It only
applies when the weight sum is 0, so all other results are unchanged.

## suave/Grid.c: grid refinement for a (practically) zero integrand

`RefineGrid` normalizes the smoothed sums of f^2 per bin with
`norm = 1/sum`. Where the integrand is (practically) zero, these sums can be
denormal (~1e-323), so `1/sum` overflows to infinity, every bin gets
`r = inf`, and the importance `((r - 1)/log(r))^1.5` is NaN. The NaN grid
then gives NaN sample coordinates.

The change keeps the previous grid in that case, as Suave already does when
the sum is exactly 0.

## suave/Fluct.c: fluctuations capped instead of pinned

`Fluct` accumulates the fluctuation of each half of a region with
`v->fluct = MaxL(f, REALL_MAX/2)`, which sets every sum to at least
`REALL_MAX/2`. All dimensions and halves then get the same fluctuation, so
Suave always bisects the widest dimension and splits the new samples evenly,
independently of the integrand; only the choice of the region to split
(largest error) is adaptive. The line is the same in the Cuba 4.2.1 and
4.2.2 releases and looks like a typo for `MinL`, an overflow cap.

The change uses `MinL`. Tested with subnucleondiffraction on 6 JIMWLK-evolved
IP-Glasma Wilson lines each for p and Pb (amplitude vs t and the
t-integrated amplitude at 1e5, 3e5, 1e6 points, against 1e7-point
references): the integrals agree within the integration error, and the
error at fixed points is about halved in most cases (no gain for the Pb
amplitude vs t, which has diffractive minima).
