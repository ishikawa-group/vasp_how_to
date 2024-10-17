# Charge analysis - Bader charge
* The atomic charge can be calculted with *Bader's method*. This charge is called **Bader charge**.
* The positive charge means atoms are positively charged (cationic), and negative charge means anionic atoms.

## INCAR
* To calculate the Bader charge, following setting in `INCAR` is needed.

* `LCHARG = .TRUE.`
* `LAECHG = .TRUE.`

* This generates the file `AECCAR0` and `AECCAR2`. These are necessary for Bader charge calculation.

---

#### その後
1. Henkelman groupのとこからbinaryかソースを落とす(http://theory.cm.utexas.edu/henkelman/code/bader/)
2. 全電子のファイルを作る
`vtstscript/chgsum.pl AECCAR0 AECCAR2` → CHGCAR_sumが出来る
3. `bader CHGCAR -ref CHGCAR_sum`とするとACF.datやらBCF.datやらが出てくる。chargeは**ACF.dat**に書いてある。valence electronの数が記述されるので、neutralだとどのくらいの電子数かをPOTCARで確認すること(ZVAL in POTCAR)
* baderはVTSTscriptの中には入っていない(=別プログラム)なので注意
* valenceのCHGCARだけでもbader analysisはできる
    * coreを入れるとNGX等を上げないと精度良くならない
    * valenceだけだとたまに異常なchargeが出るような気もする

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

#### VESTAでスピン密度をプロット
* spin-polarizedの計算をしたCHGCARからmagnetization densityを抽出
`vtstscripts/chgsplit.pl CHGCAR` --> `CHGCAR_tot`と`CHGCAR_mag`が出てくる
* `CHGCAR_mag`をVESTAで読み込みプロット <-- ファイル名をxxxCHGCARにしないと読み込めないので注意