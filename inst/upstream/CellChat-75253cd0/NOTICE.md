## Pinned CellChat source

The R delegation layer uses source files and data from `jinworks/CellChat` at commit
`75253cd0c9e68410e6e721a6d3a0419a1d7e358f` (version 2.2.0.9001). They are included so an
installed `cellchatrs` package can run without a separate CellChat checkout. The copied files
are unchanged from that commit.

CellChat is maintained by its upstream authors and is licensed under GPL-3. The upstream
license is included beside these files. The Rust inference implementation and the R shim are
the `cellchatrs` work; upstream files remain attributable to CellChat.

The following SHA-256 values pin the included files to the checkout used for this package:

```text
d7c3cadb80cd8612e5e7ba8168cb5a72c78eb6e08327ce20555477a09c628f3e  R/CellChat_class.R
575f9582b0b222e286283a45d8ed56df26295903516aef7abb44388788e800c4  R/analysis.R
fcce8c7fa5a3ae836f6184610f4e58660dd2341b62cbf8dd49301dbd5c420ac3  R/database.R
edf75f668a1527207c0d1a89316416becbbf7977a28f7226621d661366bd2cfe  R/modeling.R
061f90e679279382d96d6ccae09a825c6ad53571513df4dc4fb4ac8b9ba38895  R/utilities.R
6d80d608ad05b48fd5ec39d26099fa40cf8b389f4abb6708e27d66ab0176ffb4  R/visualization.R
582e99db2bb1b1e05d0364553893d804952bd007fef870720a9cc659a7bd3ad8  data/CellChatDB.human.rda
dd59b19646e546b3ed41f394ad2cf64cf44f3d3fcfbb58ac176f0fe1322f9413  data/CellChatDB.mouse.rda
58fff2173c93c447b3929a2f98af70d092023692c508cc0a0ec39f27505237e7  data/CellChatDB.zebrafish.rda
00fd7dc122ca83c5083892b68fd7bf006a44e029e572d6375304d2fb84e4a66c  data/PPI.human.rda
d32c5427dbddf04a0e01bd76c24feefbfb6e3c1042b7d480e0ebc8faa45f3e5e  data/PPI.mouse.rda
3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986  LICENSE
```
