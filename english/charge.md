# Charge analysis - Bader charge
* The atomic charge can be calculted with *Bader's method*. This charge is called **Bader charge**.
* The positive charge means atoms are positively charged (cationic), and negative charge means anionic atoms.

## INCAR
* It is always better to calculate the Bader charge with *all-electron charge density*, i.e. valence and core electrons.
* VASP usually doesn't give all electron charge density (valence only) so following INCAR tag is needed.
  + `LCHARG = .TRUE.`
  + `LAECHG = .TRUE.`
* This generates the files `AECCAR0` and `AECCAR2`.

## Making Bader charge sum with VTST script
### Preparation
#### VTST script
1. Download VTST scripts from Henkelman group (in University of Texas): https://theory.cm.utexas.edu/code/vtstscripts.tgz
2. Extract by `tar zxvf vtstscripts.tgz`

#### Bader
1. Download bader: https://theory.cm.utexas.edu/henkelman/code/bader/download/bader.tar.gz
2. Extract by `tar zxvf bader.tar.gz`
3. `cd bader`
4. `cp makefile.lnx_ifort Makefile`
5. `make`
6. `module load intel` (when using supercomputer)
7. Done. Command file for Bader analysis is `bader` in that directory.

### Making CHGCAR_sum
* Make all-electron CHGCAR by `vtstscript/chgsum.pl AECCAR0 AECCAR2`
  + `CHGCAR_sum` should be generated
* `bader CHGCAR -ref CHGCAR_sum`
  + Files like `ACF.dat` or `BCF.dat` are generated.
  + Bader charge is written in **ACF.dat**.
  + "CHARGE" colum gives **the number of valence electrons** for each atom in that system (so the value is "electron population" rather than "charge").
  + To calculate the charge, the number of valence electrons (in atom) should be subtracted.
  + This value is written in *ZVAL in POTCAR*.
  ```math
  {\rm charge} = -(N_{\text{electron in molecule}} - {\rm ZVAL})
  ```

#### VESTAでのCHGCARの見方
* CHGCARを開く
* "Objects" --> "Property" --> "Isosurface"で値とか色とかを調節。positiveってのが、電子をたくさん持っている(negative charge)に該当する？

#### VESTAでcharge densityを2Dプロット
* "Utilities" --> "2D Data Display"
* Sliceで平面を規定する。３つ原子を選んで書くのも良いしMiller指数でもよい
* 横のアイコンをクリックするとマウスでsliceの位置を変えれる
* 色味を調節するのはmaxとminのところで可能

#### VESTAで差電子密度
* 吸着子+表面、吸着子、表面の3つで計算しそれぞれCHGCARを用意する
* "Edit" --> "Edit Data" --> "Volumetric data" -> "Import"でsurf + ads - surfを行う -> 同様にadsも引く
* "CHGCAR"という名前でしかインポートできないのでディレクトリを別に作成する
Ref: http://renqinzhang.weebly.com/uploads/9/6/1/9/9619514/charge_density_difference.pdf
* Bader chargeのところにあるCHGCAR_sumを用いると、全電子での差電子密度もプロットできる。ただし、価電子でやるのとあまり変化はないように見える
* Mathematicaのプラグインを使った方法？(VMDでCHGをみたり、DOSを書くプラグインも説明されている)
https://www.uni-due.de/~hp0058/?file=mathplugins.html
(やってはいない)

#### Plotting the spin density with VESTA
1. Do spin-polarized calculation by setting `ISPIN = 2` and `LCHARG = .TRUE.` in `INCAR`.
2. Confirm that `CHGCAR` is generated.
3. `vtstscripts/chgsplit.pl CHGCAR` --> `CHGCAR_tot` and `CHGCAR_mag` are generated
4. Copy `CHGCAR_mag` to `tmpCHGCAR` (any xxxCHGCAR filename is OK) and load with VESTA.