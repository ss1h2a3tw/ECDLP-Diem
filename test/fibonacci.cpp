
// This test is intended to be a benchmark of computing speed

int fib(int n) {
    if(n <= 1)
        return 1;
    return fib(n-1) + fib(n-2);
}

int main() {
    return (fib(42) == 433494437) ? 0 : 1;
}
