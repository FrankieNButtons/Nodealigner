# gfa2bin-aligner

`gfa2bin-aligner` is a command line toolkit for working with variation graphs and VCF files. It extends the [`gfa2bin`](https://github.com/MoinSebi/gfa2bin) utilities with features for aligning VCF files to graph paths, extracting path information, synthesizing VCF headers and sorting variant records. The goal is to make it easy to convert coordinates from a VCF into graph path names and back again while keeping standard VCF tooling compatible.

## Features

- **Align:** Replace `#CHROM` values in a VCF using a TSV alignment of nodes to graph paths. Supports optional filtering, sorting and parallel execution.
- **Extract:** Generate a `reference.tsv` with node start/end coordinates for each path contained in a GFA file.
- **Header:** Build a valid VCF header based on keys inferred from a VCF-like file and path data from `reference.tsv`.
- **Sort:** Order VCF records by an arbitrary column (e.g. `POS`, `CHROM`, or an index).

## New in v0.0.3

- **Header subcommand:** newly added to reconstruct contig and FORMAT lines so that tools like `bcftools` accept the output without complaints.
- **Streaming alignment pipeline:** `align` can now filter unwanted chromosomes (`--skip`), normalize names (`--ignore`), optionally sort records and synthesize a header automatically unless `--no-header` is supplied.

## Installation

This repository is a Cargo project. To build the current release (version `0.0.3`) from source:

```bash
cargo build --release
```

The compiled binary will be located at `target/release/gfa2bin-aligner`.

### Typical workflow

1. **Extract** node to path coordinates from a GFA into `reference.tsv`.
2. **Align** your VCF with `reference.tsv` and a TSV alignment file to rewrite `#CHROM` to path names.
3. **Header** adds a valid VCF header using path information.
4. **Sort** re-orders the resulting VCF if required.

## Usage

Each operation is provided as a subcommand. Use `--help` on any subcommand for full argument details. Below is a quick reference to the most important flags and design choices for each mode.

### Align

```bash
gfa2bin-aligner align --vcf input.vcf --alignment alignment.tsv --reference reference.tsv
```

The alignment TSV must contain at least five tab-separated fields: `node` in the first column and `path` in the fifth column. The optional `reference.tsv` (see *Extract* below) is used as a fallback when a node is missing from the alignment file.

Key parameters:

- `--skip <str>` – comma-separated substrings. Records whose original `#CHROM` contains any substring are removed.
- `--ignore <0-5>` – normalizes `#CHROM` values. Level `0` keeps the raw name, `4` (default) restricts to `chr{1..22,X,Y,M}` and level `5` drops the `chr` prefix.
- `--sort` – sort VCF records by a column specified with `--prefix` (default: `POS`). `--reverse` reverses the order.
- `--threads <n>` – size of the Rayon thread pool. Useful for large files.
- `--no-header` – skip synthesizing a header. Without this flag `--reference` is required to create contig lines.

### Extract

```bash
gfa2bin-aligner extract --gfa graph.gfa --output reference.tsv
```

Generates a four-column `reference.tsv` (`node`, `start`, `end`, `path`) that records where each node appears in every path of the GFA. Use `--threads` to parallelize parsing on big graphs.

### Header

```bash
gfa2bin-aligner header --vcf input.vcf --reference reference.tsv --output output.vcf
```

Produces a new VCF with a valid header by merging keys inferred from the input body with path information from `reference.tsv`. Parameters mirror those of `align`:

- `--ignore <0-5>` – same normalization levels as in *Align*.
- `--threads <n>` – optional multi-threading for reading the input.
- `--output <file>` – defaults to `<input>.withheader.vcf` when omitted.

### Sort

```bash
gfa2bin-aligner sort --vcf input.vcf --prefix POS
```

Sorts a VCF by a named field or numeric index. Use `--reverse` for descending order. To keep the original file untouched, specify an explicit output name with `--output` when using `align --sort`.

## Tips

- Always run `extract` on your GFA first to obtain `reference.tsv` before aligning VCFs.
- Normalize `#CHROM` values consistently with `--ignore` to avoid mismatches between different assemblies.
- Provide a comma-separated list of contigs with `--skip` to drop unwanted chromosomes like `chrM` or scaffolds.
- Large datasets benefit from `--threads` to utilize all available CPU cores.
- After sorting, the tool inserts `.sorted` before the `.vcf` extension to prevent overwriting the unsorted output.


## TODO

- Handle gzipped VCF and TSV files directly.
- Expose the streaming combine pipeline as a dedicated subcommand.
- Publish pre-built binaries for common platforms.
- Add integration tests for typical end-to-end workflows.

