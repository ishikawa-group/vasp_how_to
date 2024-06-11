# VASPのコンパイル
* Vienna ab initio simulation package (VASP)をリモート環境でコンパイルする
* コンパイルとは、ソースコードからバイナリーファイル(実行可能ファイル)を生成すること
* 実行可能ファイルが用意されていて、それを使う場合はこの項目は飛ばしてよい

## ソースコードの入手
* VASPのソースコードを入手する
1. VASP社に料金を支払ってVASPのポータルサイトのアカウントを取得
2. ポータルサイトからvaspとpotentialの.tar.gzファイルをダウンロードしてくる

## コンパイル
1. リモート環境にvaspとpotentialのtar.gzをコピーする: `scp vasp.5.4.4.tar.gz remote:/home/your_vasp_dir`
  * remoteはリモートホスト名。`~/.ssh/config`のHostに書いてあるもの。
2. ファイルが送れたらリモート環境に移動
3. vaspを解凍: `tar zxvf vasp.5.4.4.tar.gz`
4. potentialを解凍
  ```
  mkdir -p poentials/potpaw_PBE.54
  cd potentials/potpaw_PBE.54
  mv ../../potpaw_PBE.54.tar.gz ./
  tar zxvf potpaw_PBE.54.tar.gz
  ```
4. `cd vasp.5.4.4`
5. intelコンパイラを利用可能にする: `module load intel; module load intel-mpi`
6. makefileをコピー
    * Intel one apiコンパイラ: `cp arch/makefile.include.oneapi  ./makefile.include`
    * Intelコンパイラ(古め): `cp arch/makefile.include.linux_intel  ./makefile.include`
7. コンパイル: `make all`
8. `bin`ディレクトリに`vasp_std`ができる。これがVASPの実行ファイルになる
