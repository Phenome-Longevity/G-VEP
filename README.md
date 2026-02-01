# G-VEP: GPU-Accelerated Variant Effect Prediction

GPU-accelerated variant annotation achieving 17-fold plugin speedup and 3-fold end-to-end acceleration while maintaining complete concordance with standard VEP output.

---

## Overview

Variant annotation using the Ensembl VEP plugin ecosystem is fundamentally I/O-bound, consuming over 70% of total pipeline runtime through billions of redundant disk operations. G-VEP eliminates this bottleneck by pre-computing annotation databases into sorted binary arrays and performing massively parallel binary search on GPU.

**Performance**

| Metric              | Standard VEP     | G-VEP            | Improvement     |
|:--------------------|:-----------------|:-----------------|:----------------|
| Plugin runtime      | 72 min           | 4 min            | **17× faster**  |
| End-to-end runtime  | 100 min          | 33 min           | **3× faster**   |
| Concordance         | —                | 100%             | Zero errors     |

Validated across 75 clinical whole-genome samples with zero annotation discrepancies.

---

## Installation

```bash
git clone https://github.com/emregreen/G-VEP.git
cd G-VEP
./setup.sh
```

The setup script detects existing files and downloads only what's missing.

**Requirements**

| Component          | Specification                            |
|:-------------------|:-----------------------------------------|
| Operating System   | Linux (Ubuntu 22.04+)                    |
| Docker             | Required                                 |
| Python             | 3.8+                                     |
| GPU                | NVIDIA with ≥16 GB VRAM (optional)       |
| Storage            | ~70 GB                                   |

```bash
pip install numpy pysam tqdm
pip install cupy-cuda12x  # Optional: GPU acceleration
```

Without CuPy, the pipeline falls back to CPU (slower but functional).

---

## Web Server & API

For users without local GPU infrastructure, G-VEP is available as a hosted service:

**Web Interface:** https://www.phenomeportal.org/gvep  
Upload VCF files and download annotated results through your browser.

**REST API:** Programmatic access for pipeline integration.
```bash
# Request upload URL
curl -X POST https://p7001.greenaurem.org/jobs/request-upload-url \
  -H "Content-Type: application/json" \
  -d '{"account_id": "YOUR_ID", "filename": "sample.vcf.gz"}'

# Check job status
curl "https://p7001.greenaurem.org/jobs/status/JOB_ID?account_id=YOUR_ID"
```

Full API documentation available at the web interface (API tab).

---
## Usage

```bash
./run_pipeline.sh <input.vcf.gz>
```

Output: `outputs/vep/<sample>_vep_annotated_final.vcf.gz`

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  Input VCF                                                      │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────┐
│  PASS Filter                                                     │
└──────────────────────────┬───────────────────────────────────────┘
                           │
              ┌────────────┴────────────┐
              ▼                         ▼
┌─────────────────────────┐   ┌─────────────────────────┐
│  VEP (Docker)           │   │  GPU Binary Search      │
│  · Consequence          │   │  · SpliceAI             │
│  · Transcript mapping   │   │  · dbNSFP (BayesDel)    │
│  · HGVS nomenclature    │   │  · REVEL                │
│  · LOFTEE               │   │  · AlphaMissense        │
│  · gnomAD frequencies   │   │  · dbscSNV              │
│  · ClinVar              │   │  · ClinVar (GPU)        │
└────────────┬────────────┘   └────────────┬────────────┘
             │                             │
             └──────────────┬──────────────┘
                            ▼
┌──────────────────────────────────────────────────────────────────┐
│  Merge into CSQ field                                            │
└──────────────────────────┬───────────────────────────────────────┘
                           ▼
┌──────────────────────────────────────────────────────────────────┐
│  Output VCF (standard VEP format)                                │
└──────────────────────────────────────────────────────────────────┘
```

Standard VEP plugins are I/O-bound: each plugin reads compressed databases, decompresses records, parses text, and performs string matching for every variant. G-VEP pre-computes these lookups into sorted binary indices and executes O(log n) binary search in parallel across all variants simultaneously.

---

## GPU Indices

| Database       | Annotations                                    |    Size |
|:---------------|:-----------------------------------------------|--------:|
| SpliceAI       | DS_AG, DS_AL, DS_DG, DS_DL                     |  4.5 GB |
| dbNSFP         | BayesDel_addAF, BayesDel_noAF                  |  2.0 GB |
| AlphaMissense  | Pathogenicity score, classification            |  1.1 GB |
| REVEL          | Ensemble pathogenicity score                   |  930 MB |
| dbscSNV        | ada_score, rf_score                            |  240 MB |
| ClinVar        | Clinical significance                          |    5 MB |
| **Total**      |                                                |**8.8 GB**|

The 8.8 GB index footprint fits within consumer-grade 16 GB GPUs.

---

## Variant Key Encoding

Variants are encoded as 64-bit keys using bit packing:

```
┌──────────┬──────────────┬──────────────┬──────────────┐
│ Chrom    │ Position     │ Ref hash     │ Alt hash     │
│ 8 bits   │ 28 bits      │ 14 bits      │ 14 bits      │
└──────────┴──────────────┴──────────────┴──────────────┘
```

Allele hashes use FNV-1a truncated to 14 bits. For SNVs, the four nucleotides produce distinct hash values, eliminating collision risk. Validation across 75 whole-genome samples (350M+ variant queries) showed 100% concordance.

---

## Data Components

| Component              | Description                          |      Size |
|:-----------------------|:-------------------------------------|----------:|
| GPU Indices            | Prebuilt binary search arrays        |      9 GB |
| VEP Cache              | Ensembl RefSeq cache (v113)          |     48 GB |
| LOFTEE                 | Loss-of-function plugin + data       |     13 GB |
| FASTA Reference        | GRCh38 genomic sequence              |      1 GB |
| Custom Annotations     | ClinVar, internal databases          |    300 MB |
| **Total**              |                                      | **~70 GB**|

---

## Output

The final VCF includes all standard VEP annotations with GPU-computed scores merged into the CSQ field:

| Source        | Fields                                                       |
|:--------------|:-------------------------------------------------------------|
| SpliceAI      | DS_AG, DS_AL, DS_DG, DS_DL (delta scores)                    |
| REVEL         | Score                                                        |
| dbNSFP        | BayesDel_addAF_score, BayesDel_noAF_score, predictions       |
| AlphaMissense | am_pathogenicity, am_class                                   |
| dbscSNV       | ada_score, rf_score                                          |

Output format is identical to standard VEP, ensuring compatibility with downstream clinical interpretation pipelines.

---

## Troubleshooting

**Docker permission denied**
```bash
sudo usermod -aG docker $USER
# Log out and back in
```

**CuPy installation**
```bash
nvidia-smi                    # Check CUDA version
pip install cupy-cuda11x      # CUDA 11.x
pip install cupy-cuda12x      # CUDA 12.x
```

**Verify setup**
```bash
./setup.sh                    # Re-run to check components
```

---

## Citation

If you use G-VEP in your research, please cite:

> Green E. G-VEP: GPU-Accelerated Variant Effect Prediction for Clinical Whole-Genome Sequencing. (2025)

Please also cite the underlying tools:

- **VEP:** McLaren W, et al. *Genome Biology* (2016)
- **LOFTEE:** Karczewski KJ, et al. *Nature* (2020)
- **SpliceAI:** Jaganathan K, et al. *Cell* (2019)
- **REVEL:** Ioannidis NM, et al. *Am J Hum Genet* (2016)
- **AlphaMissense:** Cheng J, et al. *Science* (2023)
- **dbNSFP:** Liu X, et al. *Genome Medicine* (2020)

---

## Author

**Emre Green**  
KTH Royal Institute of Technology & Karolinska Institute  
Phenome Longevity

---

## License

MIT