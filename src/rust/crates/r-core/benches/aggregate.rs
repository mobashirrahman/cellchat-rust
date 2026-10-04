use criterion::{criterion_group, criterion_main, Criterion};

fn placeholder(c: &mut Criterion) {
    c.bench_function("placeholder", |b| b.iter(|| std::hint::black_box(1u64)));
}

criterion_group!(benches, placeholder);
criterion_main!(benches);
