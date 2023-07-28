# Changelog

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
