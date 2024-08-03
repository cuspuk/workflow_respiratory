# Changelog

## [4.0.4](https://github.com/cuspuk/workflow_respiratory/compare/v4.0.3...v4.0.4) (2024-08-03)


### Performance Improvements

* autobump conda envs ([8f6a725](https://github.com/cuspuk/workflow_respiratory/commit/8f6a725464226838735d685b6d3fb2bf9befc653))

## [4.0.3](https://github.com/cuspuk/workflow_respiratory/compare/v4.0.2...v4.0.3) (2024-07-03)


### Performance Improvements

* autobump conda envs ([42e55f0](https://github.com/cuspuk/workflow_respiratory/commit/42e55f066b0758a77b648d5101feb1b4bbb77833))

## [4.0.2](https://github.com/cuspuk/workflow_respiratory/compare/v4.0.1...v4.0.2) (2024-06-03)


### Performance Improvements

* autobump conda envs ([1d773db](https://github.com/cuspuk/workflow_respiratory/commit/1d773db743b9b37aaf9f7508db8900dfbaa88c99))

## [4.0.1](https://github.com/cuspuk/workflow_respiratory/compare/v4.0.0...v4.0.1) (2024-05-17)


### Bug Fixes

* after repo transfer ([2eec9a4](https://github.com/cuspuk/workflow_respiratory/commit/2eec9a4a1f57bb3b687f9407ebe9e992988a437b))
* updated wrappers for new repo ([cfba283](https://github.com/cuspuk/workflow_respiratory/commit/cfba283ee6b23212c4c9a4a2bea951ee92cfd6ba))


### Performance Improvements

* bumped wrappers ([759b00d](https://github.com/cuspuk/workflow_respiratory/commit/759b00d71beacb4be9ce5250d711551361315b4a))

## [4.0.0](https://github.com/xsitarcik/respiratory/compare/v3.2.1...v4.0.0) (2024-04-16)


### ⚠ BREAKING CHANGES

* reworked to use reads module

### Features

* added depths calculation per segment ([39825dc](https://github.com/xsitarcik/respiratory/commit/39825dc72f2fc2f8a43b717fe7acc9208a99aed8))
* references are now evaluated using samtools depth ([7387f20](https://github.com/xsitarcik/respiratory/commit/7387f200ee1726c205230d4299ab7cbb6d3c8b9f))
* reworked to use reads module ([f47efc0](https://github.com/xsitarcik/respiratory/commit/f47efc02f41d446bbbe52df1094076149987fbaa))
* updated nextclade to v3 ([4b69e5c](https://github.com/xsitarcik/respiratory/commit/4b69e5c182fc2d93b74c446e74360eda68c4f14f))


### Bug Fixes

* added default production values into config ([072ab47](https://github.com/xsitarcik/respiratory/commit/072ab471f6475b980dee0b68e132133764c13bd9))
* depths json now contains only passed refs ([b662b1c](https://github.com/xsitarcik/respiratory/commit/b662b1cd9c9b6c3803dc837154034cec800b7f67))


### Performance Improvements

* bumped reads module to v3.2.1 ([0fedbcd](https://github.com/xsitarcik/respiratory/commit/0fedbcd2cb73a00657fc0a5dfb41a1c815df8d01))
* depths script marked as local ([7601488](https://github.com/xsitarcik/respiratory/commit/76014886fe467252ce450a9404d1eec543849a2e))
* removed benchmarks in rules ([b14ba4a](https://github.com/xsitarcik/respiratory/commit/b14ba4adb2a6fa4f7926904e4fb3df7aedd875d0))

## [3.2.1](https://github.com/xsitarcik/respiratory/compare/v3.2.0...v3.2.1) (2024-02-11)


### Bug Fixes

* mixed positions count no longer temporary ([71ff36d](https://github.com/xsitarcik/respiratory/commit/71ff36ddef98eb22ffa193848ad0181930aac8a5))
* turned off min mean coverage check ([0c34949](https://github.com/xsitarcik/respiratory/commit/0c349498c3e034ea20f4c6d9088dbc59e2d98066))


### Performance Improvements

* bumped wrapper versions ([337c06a](https://github.com/xsitarcik/respiratory/commit/337c06ac96e905ab29f8b983e2db97b2e16c105c))

## [3.2.0](https://github.com/xsitarcik/respiratory/compare/v3.1.7...v3.2.0) (2024-01-30)


### Features

* added another criterion for QC mapping ([e432965](https://github.com/xsitarcik/respiratory/commit/e4329653d7255812a6d205e9fed01b4770eadf47))


### Performance Improvements

* added more rules as local ([7867d70](https://github.com/xsitarcik/respiratory/commit/7867d7087a9b8c8e7590bfd57752da4e48bae584))

## [3.1.7](https://github.com/xsitarcik/respiratory/compare/v3.1.6...v3.1.7) (2023-11-29)


### Bug Fixes

* type in merged nextclades denoted by reference ([ef4bb5a](https://github.com/xsitarcik/respiratory/commit/ef4bb5a8482af633d0fc02dbbe264792f6d1306e))
* updated config to best practice values ([0cc78fb](https://github.com/xsitarcik/respiratory/commit/0cc78fbd2fbb5de96204bb7c0a9fb01cba45a5af))


### Performance Improvements

* bumped versions for tools to the latest ([51dcae2](https://github.com/xsitarcik/respiratory/commit/51dcae2338a931ff57938371890e0a02b52b8d25))

## [3.1.6](https://github.com/xsitarcik/respiratory/compare/v3.1.5...v3.1.6) (2023-11-05)


### Performance Improvements

* autobump conda envs ([f3e221f](https://github.com/xsitarcik/respiratory/commit/f3e221f82af9bbb47e761421f5b8aa0941c7e585))

## [3.1.5](https://github.com/xsitarcik/respiratory/compare/v3.1.4...v3.1.5) (2023-10-23)


### Bug Fixes

* moved reference_summary to summary dir ([de2da7c](https://github.com/xsitarcik/respiratory/commit/de2da7cbd14fd5d7cf546b230fd393b35926dfb2))

## [3.1.4](https://github.com/xsitarcik/respiratory/compare/v3.1.3...v3.1.4) (2023-10-22)


### Bug Fixes

* added full consensus for sample nextclade refs ([00185f1](https://github.com/xsitarcik/respiratory/commit/00185f17d00420dd25c6934314dc5f5bac42c5e1))

## [3.1.3](https://github.com/xsitarcik/respiratory/compare/v3.1.2...v3.1.3) (2023-10-22)


### Bug Fixes

* mapping indexes unset temporary ([035b29a](https://github.com/xsitarcik/respiratory/commit/035b29a0a53263b4f9a487e7cd5428438669e68a))

## [3.1.2](https://github.com/xsitarcik/respiratory/compare/v3.1.1...v3.1.2) (2023-10-11)


### Bug Fixes

* added retry for fallible krona update ([58c4152](https://github.com/xsitarcik/respiratory/commit/58c4152309597b01b68ae1781bbcf154bf97142c))
* passed bams are copied with indexes ([95d1b0a](https://github.com/xsitarcik/respiratory/commit/95d1b0a8f4deda61b2d13b40156f89a3119a1c7b))
* qual output in ivar consensus set to temporary ([ce052da](https://github.com/xsitarcik/respiratory/commit/ce052da659de464f0a1d5b265b4e86aefce22129))
* reformatted mixed positions counts summary file ([ed414eb](https://github.com/xsitarcik/respiratory/commit/ed414eb8e0ad0ad47b24b1c5f7cfbf07ac35030b))
* uncompressed files during decontamination set to temporary ([59cf7f4](https://github.com/xsitarcik/respiratory/commit/59cf7f42b87a43c20e698b7a36d66760804913dd))

## [3.1.1](https://github.com/xsitarcik/respiratory/compare/v3.1.0...v3.1.1) (2023-10-10)


### Bug Fixes

* separated passed bams from all bams ([6f6bcb0](https://github.com/xsitarcik/respiratory/commit/6f6bcb0e35db0c50da102520b5de138571197e14))

## [3.1.0](https://github.com/xsitarcik/respiratory/compare/v3.0.2...v3.1.0) (2023-10-08)


### Features

* added conversion of nextclade tsv into html ([b2985b8](https://github.com/xsitarcik/respiratory/commit/b2985b872cb69a14425731a6257ec859ba80cd8d))
* aggregated sample consensuses for each reference ([eae9f56](https://github.com/xsitarcik/respiratory/commit/eae9f560c3551626388e8c10a87dc3df64752c58))
* merged all nextclade tsvs into one ([0c55c73](https://github.com/xsitarcik/respiratory/commit/0c55c733bc10b46e2230172822debe6331c4bb6a))
* merging nextclade tsv per sample ([3b7870d](https://github.com/xsitarcik/respiratory/commit/3b7870d9c549526a3ae44aafe8c37fe80809717d))


### Bug Fixes

* aggregation runs when more than one sample ([9be8079](https://github.com/xsitarcik/respiratory/commit/9be807961f23ac99bd6e519df81ecf70ea7fa0db))
* ivar variants produced also in vcf format ([b178dbe](https://github.com/xsitarcik/respiratory/commit/b178dbefe9636bf13d0a388df612c48a25ee0a88))
* moved nextclade outputs into own subdirectory ([657293a](https://github.com/xsitarcik/respiratory/commit/657293a8832d50f517f0f4fad51b1a784453817e))
* reports grouped per sample ([5848089](https://github.com/xsitarcik/respiratory/commit/58480899caeed7f32349cf876f593dd482ecd621))

## [3.0.2](https://github.com/xsitarcik/respiratory/compare/v3.0.1...v3.0.2) (2023-10-03)


### Bug Fixes

* fixed validation schema for consensus threshold ([e64a176](https://github.com/xsitarcik/respiratory/commit/e64a176727800f67f134e8235c9a2b383c027097))

## [3.0.1](https://github.com/xsitarcik/respiratory/compare/v3.0.0...v3.0.1) (2023-08-16)


### Bug Fixes

* fixing parsing paired cutadapt params ([3d38e39](https://github.com/xsitarcik/respiratory/commit/3d38e39c47ce7e390497915cf64fa3cfdb4934db))

## [3.0.0](https://github.com/xsitarcik/respiratory/compare/v2.0.2...v3.0.0) (2023-08-14)


### ⚠ BREAKING CHANGES

* input kraken path required to be already tagged

### Features

* added nextclade using accession and version tag ([9ab76c7](https://github.com/xsitarcik/respiratory/commit/9ab76c739c9e93b00620dc944d3c49005b24ee51))
* input kraken path required to be already tagged ([9f525af](https://github.com/xsitarcik/respiratory/commit/9f525afebbe333c94650c6c6900e02a9dc720933))
* kraken input merged into path with tag ([fe03f10](https://github.com/xsitarcik/respiratory/commit/fe03f10b49c2da4d671feda5b12d4a285ad847e0))


### Bug Fixes

* nextclade are downloaded into common directory ([d4085e6](https://github.com/xsitarcik/respiratory/commit/d4085e6827242ed95a4ca39f62afa934e3952456))

## [2.0.2](https://github.com/xsitarcik/respiratory/compare/v2.0.1...v2.0.2) (2023-07-28)


### Bug Fixes

* added configurable kraken mode to save memory ([e2b446a](https://github.com/xsitarcik/respiratory/commit/e2b446a80e4e1ffc8c22e1e208d2f19bc14bce89))

## [2.0.1](https://github.com/xsitarcik/respiratory/compare/v2.0.0...v2.0.1) (2023-07-28)


### Bug Fixes

* added configurable krona directory in the config ([37fd19b](https://github.com/xsitarcik/respiratory/commit/37fd19b1c513fa082da067f51aa088b363a17e3a))
* removed decontamination rule stdout in log ([77f604d](https://github.com/xsitarcik/respiratory/commit/77f604d352fa416c7b302909548d2eceef1db4e6))

## [2.0.0](https://github.com/xsitarcik/respiratory/compare/v1.2.0...v2.0.0) (2023-07-26)


### ⚠ BREAKING CHANGES

* changed data input to use pep files

### Features

* changed data input to use pep files ([1be145c](https://github.com/xsitarcik/respiratory/commit/1be145cc2d1278b9a8d714d6e65496fda1f2c876))

## [1.2.0](https://github.com/xsitarcik/respiratory/compare/v1.1.3...v1.2.0) (2023-07-17)


### Features

* nextclade supports multiple segments of reference ([7ad5608](https://github.com/xsitarcik/respiratory/commit/7ad56083d3455ffc94954af267b86f018a90ee95))


### Bug Fixes

* nextclade uses default tag if not provided ([6ce0a28](https://github.com/xsitarcik/respiratory/commit/6ce0a281663fc623c48f56d5e864bcaa64eeb1af))

## [1.1.3](https://github.com/xsitarcik/respiratory/compare/v1.1.2...v1.1.3) (2023-07-04)


### Bug Fixes

* fixed parsing refernce name from qualimap dir ([d895c13](https://github.com/xsitarcik/respiratory/commit/d895c13f40a9bec35e4294d70d3dcc839debd134))

## [1.1.2](https://github.com/xsitarcik/respiratory/compare/v1.1.1...v1.1.2) (2023-07-04)


### Bug Fixes

* change bamqc to different wrapper ([abb2b86](https://github.com/xsitarcik/respiratory/commit/abb2b86ad42603ac3bd763bc2a48293532e3273a))

## [1.1.1](https://github.com/xsitarcik/respiratory/compare/v1.1.0...v1.1.1) (2023-07-03)


### Bug Fixes

* fixed parsing reference name in nonempty bams ([9a7f584](https://github.com/xsitarcik/respiratory/commit/9a7f58436e8a55ee6561fa72e9d905b80d8dbda4))

## [1.1.0](https://github.com/xsitarcik/respiratory/compare/v1.0.5...v1.1.0) (2023-07-03)


### Features

* added resource handling from config ([8749ba6](https://github.com/xsitarcik/respiratory/commit/8749ba6ce02ed35e0de76d3fc543869ed1740cc0))

## [1.0.5](https://github.com/xsitarcik/respiratory/compare/v1.0.4...v1.0.5) (2023-06-30)


### Bug Fixes

* picard changed wrapper ([d0e64af](https://github.com/xsitarcik/respiratory/commit/d0e64af476bd90078700063c2e4dc7c2a398bf84))

## [1.0.4](https://github.com/xsitarcik/respiratory/compare/v1.0.3...v1.0.4) (2023-06-30)


### Bug Fixes

* added biopython into krakentools env ([6ce5f6b](https://github.com/xsitarcik/respiratory/commit/6ce5f6b60a94bee241fec2bf2c0c13ad8c48fdf3))

## [1.0.3](https://github.com/xsitarcik/respiratory/compare/v1.0.2...v1.0.3) (2023-06-30)


### Bug Fixes

* added cutadapt to replace trimmomatic ([1e56a2a](https://github.com/xsitarcik/respiratory/commit/1e56a2acc2d0726c4d86b741cb145252233814bc))

## [1.0.2](https://github.com/xsitarcik/respiratory/compare/v1.0.1...v1.0.2) (2023-06-29)


### Bug Fixes

* added memory resource handling for trimming ([455a78a](https://github.com/xsitarcik/respiratory/commit/455a78ad109939555dd67bca934ac22ffba76b06))

## [1.0.1](https://github.com/xsitarcik/respiratory/compare/v1.0.0...v1.0.1) (2023-06-26)


### Bug Fixes

* added config parameter for metadata to references ([d00a368](https://github.com/xsitarcik/respiratory/commit/d00a368c14e627deefd41a3cc78be4b32b09a2ed))
* added parameter for kraken dir to control downloading kraken database ([91f105d](https://github.com/xsitarcik/respiratory/commit/91f105d7c5faa00cc8e094d1f6bbfd237b9b2e3f))
* added sample_dirpath parameter to allow analyse sample from anywhere ([030f433](https://github.com/xsitarcik/respiratory/commit/030f433c9663a5bc3d3385116cca95ed50168015))
* mapping results reorganized to have sample as the upper level ([32e360c](https://github.com/xsitarcik/respiratory/commit/32e360ca1d65afa214ec963a7da138a5576ba55a))

## 1.0.0 (2023-06-19)


### Features

* Initial version of functioning workflow ([#1](https://github.com/xsitarcik/respiratory/issues/1)) ([b0893d4](https://github.com/xsitarcik/respiratory/commit/b0893d4b83960973015de4360dfd71fea491a909))
