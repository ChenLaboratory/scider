# Hash two 15-bytes signed integers into one 32-bytes integer.

Hash two 15-bytes signed integers into one 32-bytes integer.

## Usage

``` r
coord_hash(a, b)
```

## Arguments

- a, b:

  Integer vectors of same lengths. Can be negative.

## Details

Should work for a,b in range (-2^14, 2^14-1) which is good enough for
our purpose.

Integer in R is 32-bytes but reserve 1 byte for NA
