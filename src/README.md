# インプットファイルの説明（EnKF.inp）

インプットファイル4DVar.inpの中身は以下のようになっています．
```
2                  !--- Mode (1:TLM&ADJ Check, 2:4DV, 3:4DV restart)
1                  !--- Problem 1:Karman vortex, 2:Vortex advection
17 50              !--- Number of 4DV cycles, iteration on each cycle
0.0                !--- Error variance of pseudo measurement
0.1                !--- Measurement error variance
0.0                !--- Coefficient of background term
4 4 40             !--- Every iskip & jskip & tskip for measuremnt 
-2.6d0 6.0d0       !--- Left and right of measurement area
-2.6d0 2.6d0       !--- Bottom and top of measurement area
------------------------------------------------------------------------------
150 150            !--- DA window step, Total number of time step
100                !--- Reynolds number 
80                 !--- Number of mesh: jmax (imax x jmax, imax = 3*jmax)
0.02               !--- time step
1                  !--- 0:initial, 1:restart
1                  !--- 0:No output to screen,1:otherwise
10                 !--- Output plot3d interval
20                 !--- Output history interval
```

# ソースコードの説明

## EnKF_main.f90

メインルーチンです．インプットファイル（EnKF.inp）の読み込みから，アンサンブルカルマンフィルタの処理まで全てを行っています．

## m_random3.f90

平均0，分散1の正規乱数ベクトルを出力するサブルーチンです．変数の数n_varとアンサンブル数n_ensを入力して，そのサイズの正規乱数ベクトルを得ます．

## m_ranmean3.f90

変数の数n_varとアンサンブル数n_ensの配列を入力して，アンサンブル平均を計算します．

## m_ranvar3.f90

変数の数n_varとアンサンブル数n_ensの配列およびアンサンブル平均ベクトルを入力して，分散を計算します．

## makefile

ソースコードをコンパイル・リンクして，実行ファイル4dvarを作成します．通常はsrcフォルダ内で`make`と入力します．デバッグモードでコンパイルする場合には`make debug`とします．

## mod_variables.f90

使用する変数を定義しています．

## sub_kfilter.f90

アンサンブルカルマンフィルタによるフィルタリングを行うサブルーチンです．アンサンブルメンバーの含まれる配列と観測ベクトルの配列を入力すると，アンサンブルが更新されて出力されます．カルマンゲインにおける逆行列計算に関して，`n_mes x n_mes`行列を扱うサブルーチンと，`n_ens x n_ens`で計算するサブルーチンが含まれています．これらの切り替えはインプットファイルEnKF.inpで行います．

## sub_measure.f90

観測モデルを定義しています．どの範囲の計算格子点を何点おきに計測点とするかを定義してます．

## sub_bc_outer.f90

計算領域の外側境界に流れの流入・流出などの境界条件を与えます．

## sub_bc_wall.f90

無限長さ角柱（2次元計算領域における正方形）まわりの境界条件を埋め込み境界法によって与えるサブルーチンです．

## sub_hsmac.f90

圧力および速度の同時過緩和をHSMAC法によって行います．

## sub_initial.f90

使用する配列の確保，計算領域の定義を行います．sub_addvtx(xc,yc)は計算領域内にBurnham-Hallock渦を1つ設置するサブルーチンです．リスタートファイルの入力を行うサブルーチンも含まれています．

## sub_plot3d.f90

流れ場の可視化用にPlot3dファイル形式で結果を出力します．Plot3dファイル（計算格子：mesh.g, 結果ファイル：*.q）はTecplotやParaviewで可視化できます．

## sub_rhs3rd.f90

流れ場の時間発展をオイラー陽解法（1次精度陽解法）で行います．対流項はKuwahara-Kawamuraスキーム，粘性項は2次精度中心差分で離散化します．

## sub_utils.f90

アンサンブルを生成するためのサブルーチンが含まれています．


# 初回コンパイル時にダウンロードされるソースコード

## svd.f

m_kfilter.f90内で逆行列を計算する際に利用します．本リポジトリには含まれませんが，初回の`make`実行時に以下サイトから自動的にダウンロードされます．

EISPACKライブラリ  
http://www.netlib.org/eispack/svd.f

## pythag.f

m_kfilter.f90内で逆行列を計算する際に利用します．本リポジトリには含まれませんが，初回の`make`実行時に以下サイトから自動的にダウンロードされます．

EISPACKライブラリ  
http://www.netlib.no/netlib/eispack/3090vf/double/pythag.f
