#

## Legend

```mermaid
flowchart LR;
  subgraph LEG[LEGEND];
  L1[Output];
  L2[Temporary Output];

  style L1 fill:#0f0;
  style L2 fill:#ff0;
  end;


```

## Diagram

```mermaid
flowchart TB;
  classDef OUTPUT fill:#0f0;
  classDef TEMP fill:#ff0;

  A[/Raw reads/];
  B[Trimmed reads]:::TEMP;
  F[Decontaminated reads]:::TEMP;
  C[Kraken mapping result]:::TEMP;
  E[Krona report]:::OUTPUT;
  D[FastQC report]:::OUTPUT;
  Z[Summary]:::OUTPUT;
  H[Mapping quality summary]:::OUTPUT;

  style S1_1 fill:#0f0;
  style S2_0 fill:#0f0;
  style S3_2 fill:#0f0;
  style S3_0 fill:#0f0;
  style S3_1 fill:#0f0;
  style S3_3 fill:#0f0;

  A-->|trimmomatic|B;
  B-->|kraken|C;
  C-->E;
  B-->F;
  C-->F;
  F-->|bwa|S1;
  B-->|FastQC|D;

  S1-->|bamQC|S2;
  S2-->|summarize qualimap|H;
  H--->|evaluate coverage|S3;
  S1-->|ivar consensus|S3_S1;
  S1-->|ivar variants|S3_1;
  S3--->|Summarize|Z;

  O>config params]-.->S3_2;

  subgraph S1[For every reference sequence in the panel];
  direction TB;
    S1_0[Mapping output]:::TEMP;
    S1_1[Deduplicated BAM];
    S1_0-->|picard markDuplicates|S1_1;
  end;

  subgraph S2[For every non-empty BAM]
  direction TB;
    S2_0[Qualimap report];
  end;

  subgraph S3[For every BAM passing quality check]
  direction TB;
    S3_D{reference};
    S3_0[Consensus];
    S3_1[Variants tsv];
    S3_2[Mixed positions];
    S3_3[Nextclade report];
    S3_S1-->|concatenate|S3_0;
    S3_1-->|count|S3_2;
    S3_D-->|is in nextclade DB|S3_3;
    S3_0-->|nextclade|S3_3;
  end;

  subgraph S3_S1[For every segment of reference];
  direction TB;
    S3_S1_0[Segment's consensus]:::TEMP;
  end;
```
