use criterion::{black_box, criterion_group, criterion_main, Criterion};
use splitter::simdstuffs::mismatch_count;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("mismatch_comp", |b| {
        b.iter(|| mismatch_count(black_box(&[0; 32]), black_box(&[0; 32]), 0))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
