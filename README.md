# Implementing Diem's variation of Index Calculus in C++

Paper: https://eprint.iacr.org/2015/1022.pdf

This algorithm is to solve discrete log problem in elliptic curve with Index Calculus.

## Building and Testing

### Step 1
```shell
mkdir build
cd build
```

### Step 2
```shell
cmake ..
```
Or if you want to build in release mode:
```shell
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Step 3
```shell
make
```

### Step 4
```shell
make test
```

