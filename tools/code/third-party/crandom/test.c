int main(void) {
    double a;
    init_crandom(10);
    get_crandom(&a);
    printf("extracted random %lf ",a);
    return 0;
};

