##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
##  Factor analysis with lassoにおけるGtkオブジェクト
##
##  ファイル名：plot.fanc.R
##  ファイル内容：
##  作成者：YAMAMOTO, Michio
##  作成日：2012年02月02日
##  最終更新日：2012年03月19日
##  コメント："X", "F"の表示を削除（031012）
##           とりあえず使える形にした（031012）
##           マイナスの負荷量は赤で表すこととした（031012）
##           gammaのバーを追加，因子数が奇数の場合のプロット修正（031912）
##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★

##-----コールバック関数定義-----##
##expose-eventシグナルに対するコールバック関数
##図形描画関数
cbExposeCanvas <- function (gtk.widget, data)
{
  N.var <- info.fanc$N.var
  N.fac <- info.fanc$N.fac
  Window.Width <- info.fanc$Window.Width
  Window.Height <- info.fanc$Window.Height
  Rad.Ellipse <- info.fanc$Rad.Ellipse
  Len.Rec <- info.fanc$Len.Rec

  drawable <- gtk.widget$window
  HN.var <- N.var / 2
  HN.fac <- N.fac / 2

  cr <- gdkCairoCreate (drawable)
  cr.t <- gdkCairoCreate (drawable)
  Mat <- cairoMatrixInit (0, 0, 0, 0, 0, 0)$matrix

  ##-------------------
  ##   因子の楕円を描画
  ##-------------------
  Rem.N.fac <- N.fac %% 2
  cairoTranslate (cr, Window.Width / 2, Rad.Ellipse + 20) ##原点の移動
  cairoTranslate (cr.t, Window.Width / 2, Rad.Ellipse + 20) ##原点の移動

  cairoSetLineWidth (cr, 2.5)

  cairoScale (cr, 1.0, 0.5) ##楕円のスケール

  if (Rem.N.fac == 0) { ##因子数が偶数の場合
    for (i in 1:HN.fac) {
      ii <- i - 1

      ##楕円の描画
      cairoArc (cr, -Rad.Ellipse * 5 / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 0,
                Rad.Ellipse, 0.0, 2.0 * pi)
      cairoStroke (cr)
      cairoArc (cr,  Rad.Ellipse * 5 / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 0,
                Rad.Ellipse, 0.0, 2.0 * pi)
      cairoStroke (cr)

      cairoSetFontSize (cr.t, 20)

      ##ラベルの描画
      text <- sprintf ("f%d", HN.fac - ii)
      cairoMoveTo (cr.t, -Rad.Ellipse * 5 / 4 - Rad.Ellipse / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 5)
      cairoShowText (cr.t, text)

      text <- sprintf ("f%d", HN.fac + 1 + ii)
      cairoMoveTo (cr.t,  Rad.Ellipse * 5 / 4 - Rad.Ellipse / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 5)
      cairoShowText (cr.t, text)
    }
  }
  else if (Rem.N.fac != 0) { ##因子数が奇数の場合
    cairoArc (cr, 0.0, 0.0, Rad.Ellipse, 0.0, 2.0 * pi)
    cairoStroke (cr)
    cairoSetFontSize (cr.t, 20)
    text <- sprintf ("f%d", floor(HN.fac) + 1)
    cairoMoveTo (cr.t, 0 - Rad.Ellipse / 4, 5)
    cairoShowText (cr.t, text)

    if (floor(HN.fac) != 0) {
      for (i in 1:floor(HN.fac)) {
        ii <- i - 1

        ##楕円の描画
        cairoArc (cr, -Rad.Ellipse * 10 / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 0,
                  Rad.Ellipse, 0.0, 2.0 * pi)
        cairoStroke (cr)
        cairoArc (cr,  Rad.Ellipse * 10 / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 0,
                  Rad.Ellipse, 0.0, 2.0 * pi)
        cairoStroke (cr)

        cairoSetFontSize (cr.t, 20)

        ##ラベルの描画
        text <- sprintf ("f%d", floor(HN.fac) - ii)
        cairoMoveTo (cr.t, -Rad.Ellipse * 10 / 4 - Rad.Ellipse / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 5)
        cairoShowText (cr.t, text)

        text <- sprintf ("f%d", floor(HN.fac) + ii + 2)
        cairoMoveTo (cr.t,  Rad.Ellipse * 10 / 4 - Rad.Ellipse / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 5)
        cairoShowText (cr.t, text)
      }
    }
  } ##因子の楕円の描画終了

  ##これまでの座標系の変換を元に戻す
  cairoGetMatrix (cr, Mat)
  cairoMatrixInvert (Mat)
  cairoTransform (cr, Mat)
  cairoGetMatrix (cr.t, Mat)
  cairoMatrixInvert (Mat)
  cairoTransform (cr.t, Mat)


  ##-------------------
  ##  変数の四角を描画
  ##-------------------
  ##変数の数の偶奇で場合わけ
  Rem.N.var <- N.var %% 2
  cairoTranslate (  cr, Window.Width / 2, Window.Height - Len.Rec * 5 / 2)
  cairoTranslate (cr.t, Window.Width / 2, Window.Height - Len.Rec * 5 / 2)

  cairoSetLineWidth (cr, 2.5)

  if (Rem.N.var == 0) {##変数の数が偶数の場合
    for (i in 1:HN.var) {
      ii <- i - 1
      cairoRectangle (cr, -Len.Rec * 5 / 4 - ii * Len.Rec * 3 / 2, 0, Len.Rec, Len.Rec)
      cairoStroke (cr)
      cairoRectangle (cr, Len.Rec / 4 + ii * Len.Rec * 3 / 2, 0, Len.Rec, Len.Rec)
      cairoStroke (cr)

      cairoSetFontSize (cr.t, 20)

      text <- sprintf ("x%d", floor(HN.var) - ii)
      cairoMoveTo (cr.t, -Len.Rec * 5 / 4 - ii * Len.Rec * 3 / 2 + Len.Rec / 4, Len.Rec * 2 / 3)
      cairoShowText (cr.t, text)

      text <- sprintf ("x%d", floor(HN.var) + 1 + ii)
      cairoMoveTo (cr.t, Len.Rec / 4 + ii * Len.Rec * 3 / 2 + Len.Rec / 4, Len.Rec * 2 / 3)
      cairoShowText (cr.t, text)
    }
  }
  else if (Rem.N.var != 0) { ##変数の数が奇数の場合
    cairoRectangle (cr, -Len.Rec / 2, 0, Len.Rec, Len.Rec)
    cairoStroke (cr)
    cairoSetFontSize (cr.t, 20)
    text <- sprintf("x%d", floor(HN.var) + 1)
    cairoMoveTo (cr.t, 0 - Len.Rec / 2 + Len.Rec / 3, Len.Rec * 2 / 3)
    cairoShowText (cr.t, text)

    for (i in 1:HN.var) {
      ii <- i - 1
      cairoRectangle (cr, -Len.Rec * 2 - ii * Len.Rec * 3 / 2, 0, Len.Rec, Len.Rec)
      cairoStroke (cr)
      cairoRectangle (cr, Len.Rec      + ii * Len.Rec * 3 / 2, 0, Len.Rec, Len.Rec)
      cairoStroke (cr)

      cairoSetFontSize (cr.t, 20)

      text <- sprintf ("x%d", floor(HN.var) - ii)
      cairoMoveTo (cr.t, -Len.Rec * 2 - ii * Len.Rec * 3 / 2 + Len.Rec / 3 - 0.5, Len.Rec * 2 / 3)
      cairoShowText (cr.t, text)

      text <- sprintf ("x%d", floor(HN.var) + 2 + ii)
      cairoMoveTo (cr.t, Len.Rec      + ii * Len.Rec * 3 / 2 + Len.Rec / 3 - 0.5, Len.Rec * 2 / 3)
      cairoShowText (cr.t, text)
    }
  }

  ##これまでの座標系の変換を元に戻す
  cairoGetMatrix (cr, Mat)
  cairoMatrixInvert (Mat)
  cairoTransform (cr, Mat)
  cairoGetMatrix (cr.t, Mat)
  cairoMatrixInvert (Mat)
  cairoTransform (cr.t, Mat)



  ##---------------------------
  ##  因子と変数を結ぶ線を描画
  ##---------------------------
  if (Rem.N.fac == 0) { ##因子数が偶数の場合
    cr.Fac.X <- -Rad.Ellipse * 5 / 4 - (HN.fac - 1) * 2 * Rad.Ellipse * 5 / 4 + Window.Width / 2
    cr.Fac.Y <- Rad.Ellipse + 20 + Rad.Ellipse / 2
  }
  else if (Rem.N.fac != 0) { ##因子数が奇数の場合
    cr.Fac.X <- -Rad.Ellipse * 10 / 4 - (floor(HN.fac) - 1) * 2 * Rad.Ellipse * 5 / 4 + Window.Width / 2
    cr.Fac.Y <- Rad.Ellipse + 20 + Rad.Ellipse / 2
  }

  if (Rem.N.var == 0) { ##変数の数が偶数の場合
    cr.Var.X <- -Len.Rec * 5 / 4 - (HN.var - 1) * Len.Rec * 3 / 2 + Window.Width / 2 + Len.Rec / 2
    cr.Var.Y <- Window.Height - Len.Rec * 5 / 2
  }
  else if (Rem.N.var != 0) { ##変数の数が奇数の場合
    cr.Var.X <- -Len.Rec * 2 - (floor(HN.var) - 1) * Len.Rec * 3 / 2 + Window.Width / 2 + Len.Rec / 2
    cr.Var.Y <- Window.Height - Len.Rec * 5 / 2
  }

  ##全てのパスを個々に描画する
  cairoSetLineJoin (cr, 'CAIRO_LINE_JOIN_BEVEL')
##  cairoSetLineJoin (cr, 'CAIRO_LINE_JOIN_MITER')
  for (i in 1:N.var) {
    ii <- i - 1
    for (j in 1:N.fac) {
      jj <- j - 1
      cairoSetLineWidth (cr, abs(info.fanc$L[i, j, info.fanc$num.lambda, info.fanc$num.gamma] * 10))

      if (info.fanc$L[i, j, info.fanc$num.lambda, info.fanc$num.gamma] < 0) {
        cairoSetSourceRgb (cr, 255.0, 0.0, 0.0)
      }
      else {
        cairoSetSourceRgb (cr, 0.0, 0.0, 0.0)
      }
      cairoMoveTo (cr, cr.Fac.X + jj * 2 * Rad.Ellipse * 5 / 4, cr.Fac.Y)
      cairoLineTo (cr, cr.Var.X + ii * Len.Rec * 3 / 2, cr.Var.Y)
      cairoStroke (cr)
    }
  }
}



##GtkWidget labelのexpose-eventシグナルに対するコールバック関数
##-----for 'rho'-----
cbExposeLabelLambda <- function (gtk.widget, data)
{
  drawable <- gtk.widget$window

  cr <- gdkCairoCreate (drawable)

  text <- sprintf ("rho : %f", info.fanc$lambda.current)
  cairoSetFontSize (cr, 15)
  cairoMoveTo (cr, 0, 30)
  cairoShowText (cr, text)
}

##-----for 'GFI'-----
cbExposeLabelGFI <- function (gtk.widget, data)
{
  drawable <- gtk.widget$window

  cr <- gdkCairoCreate (drawable)

  text <- sprintf ("GFI : %f", info.fanc$GFI.current)
  cairoSetFontSize (cr, 15)
  cairoMoveTo (cr, 0, 30)
  cairoShowText (cr, text)
}

##-----for 'AGFI'-----
cbExposeLabelAGFI <- function (gtk.widget, data)
{
  drawable <- gtk.widget$window

  cr <- gdkCairoCreate (drawable)

  text <- sprintf ("AGFI : %f", info.fanc$AGFI.current)
  cairoSetFontSize (cr, 15)
  cairoMoveTo (cr, 0, 30)
  cairoShowText (cr, text)
}

##-----for 'AIC'-----
cbExposeLabelAIC <- function (gtk.widget, data)
{
  drawable <- gtk.widget$window

  cr <- gdkCairoCreate (drawable)

  text <- sprintf ("AIC : %f", info.fanc$AIC.current)
  cairoSetFontSize (cr, 15)
  cairoMoveTo (cr, 0, 30)
  cairoShowText (cr, text)
}

##-----for 'BIC'-----
cbExposeLabelBIC <- function (gtk.widget, data)
{
  drawable <- gtk.widget$window

  cr <- gdkCairoCreate (drawable)

  text <- sprintf ("BIC : %f", info.fanc$BIC.current)
  cairoSetFontSize (cr, 15)
  cairoMoveTo (cr, 0, 30)
  cairoShowText (cr, text)
}

##-----for 'CAIC'-----
cbExposeLabelCAIC <- function (gtk.widget, data)
{
  drawable <- gtk.widget$window

  cr <- gdkCairoCreate (drawable)

  text <- sprintf ("CAIC : %f", info.fanc$CAIC.current)
  cairoSetFontSize (cr, 15)
  cairoMoveTo (cr, 0, 30)
  cairoShowText (cr, text)
}


##-----for 'gamma'-----
cbExposeLabelGamma <- function (gtk.widget, data)
{
  drawable <- gtk.widget$window

  cr <- gdkCairoCreate (drawable)

  text <- sprintf ("gam : %f", info.fanc$gamma.current)
  cairoSetFontSize (cr, 15)
  cairoMoveTo (cr, 0, 30)
  cairoShowText (cr, text)
}


##GtkScaleウィジェット（"canvas"）のvalue-changedシグナルに対するコールバック関数
##-----for 'rho'-----
cbValueChangedLambda <- function (gtk.scale, data)
{
  if (gtkRangeGetValue (gtk.scale) == 0) {
    info.fanc$num.lambda <<- 1
  }
  else {
    info.fanc$num.lambda <<- ceiling(gtkRangeGetValue (gtk.scale) * info.fanc$N.lambda) ##誤差を修正
  }
  gtkWidgetQueueDraw (data) ##canvasに対してexpose-eventシグナルを出す
}

##-----for 'gamma'-----
cbValueChangedGamma <- function (gtk.scale, data)
{
  if (gtkRangeGetValue (gtk.scale) == 0) {
    info.fanc$num.gamma <<- 1
  }
  else {
    info.fanc$num.gamma <<- ceiling(gtkRangeGetValue (gtk.scale) * info.fanc$N.gamma) ##誤差を修正
  }

  info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]

  gtkWidgetQueueDraw (data) ##canvasに対してexpose-eventシグナルを出す
}



##GtkScaleウィジェット（"label"）のvalue-changedシグナルに対するコールバック関数
##-----for 'rho'-----
cbValueChangedLabelLambda <- function (gtk.scale, data)
{
  #if (gtkRangeGetValue (gtk.scale) == 0) {
  #  info.fanc$lambda.current <<- info.fanc$lambdas[1, info.fanc$num.gamma]
  #}
  #else {
  #  info.fanc$lambda.current <<- info.fanc$lambdas[ceiling(gtkRangeGetValue (gtk.scale) * info.fanc$N.lambda), info.fanc$num.gamma] ##誤差を修正
  #}
  info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$GFI.current <<- info.fanc$GFIs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$AGFI.current <<- info.fanc$AGFIs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$AIC.current <<- info.fanc$AICs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$BIC.current <<- info.fanc$BICs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$CAIC.current <<- info.fanc$CAICs[info.fanc$num.lambda, info.fanc$num.gamma]

  gtkWidgetQueueDraw (data[[1]]) ##label.lambdaに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[2]]) ##label.GFIに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[3]]) ##label.AGFIに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[4]]) ##label.AICに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[5]]) ##label.BICに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[6]]) ##label.CAICに対してexpose-eventシグナルを出す
}

##-----for 'gamma'-----
cbValueChangedLabelGamma <- function (gtk.scale, data)
{
  if (gtkRangeGetValue (gtk.scale) == 0) {
    info.fanc$gamma.current <<- info.fanc$gammas[1]
  }
  else {
    info.fanc$gamma.current <<- info.fanc$gammas[ceiling(gtkRangeGetValue (gtk.scale) * info.fanc$N.gamma)] ##誤差を修正
  }

  info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$GFI.current <<- info.fanc$GFIs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$AGFI.current <<- info.fanc$AGFIs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$AIC.current <<- info.fanc$AICs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$BIC.current <<- info.fanc$BICs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$CAIC.current <<- info.fanc$CAICs[info.fanc$num.lambda, info.fanc$num.gamma]

  gtkWidgetQueueDraw (data[[1]]) ##label.gammaに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[2]]) ##label.lambdaに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[3]]) ##label.GFIに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[4]]) ##label.AGFIに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[5]]) ##label.AICに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[6]]) ##label.BICに対してexpose-eventシグナルを出す
  gtkWidgetQueueDraw (data[[7]]) ##label.CAICに対してexpose-eventシグナルを出す
}



##-----描画するためのウィンドウ本体-----##
MakeInterface <- function (gchar.title)
{
  ##"info" includes 'L' and 'lambdas'.
  ##'L' is a p*m*r array of pattern matrices
  ##'lambdas' is a r*1 vector of tuning parameter

  lambdas <- info.fanc$lambdas
  gammas <- info.fanc$gammas

  ##GTK+上ではlambda, gammaの範囲を(0, 1)に限定する
  ##for 'rho'
  N.lambda <- length(lambdas[,1])
  Min.lambda <- 0
  Max.lambda <- 1
  Step.lambda <- 1 / (N.lambda - 1)

  ##for 'gamma'
  N.gamma <- length(gammas)
  Min.gamma <- 0
  Max.gamma <- 1
  Step.gamma <- 1 / (N.gamma - 1)

  window <- gtkWindowNew (show=TRUE)
  gtkWindowSetTitle (window, gchar.title)
##  gtkWidgetSetSizeRequest (window, info.fanc$Window.Width, info.fanc$Window.Height)
  gtkWidgetSetSizeRequest (window, info.fanc$Window.Width, -1)
  gtkContainerSetBorderWidth (window, 5)
##  gSignalConnect (window, "destroy", gtkMainQuit, NULL) ##なくても良いかもしれない


  ##-------------------
  ##    パス図の描画
  ##-------------------
  vbox <- gtkVBoxNew (FALSE, 3)
  gtkContainerAdd (window, vbox)
  gtkContainerSetBorderWidth (vbox, 5)

  alignment <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (vbox, alignment, TRUE, TRUE, 0)

  canvas <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment, canvas)
  gtkWidgetSetSizeRequest (canvas, info.fanc$Window.Width, info.fanc$Window.Height)
  gSignalConnect (canvas, "expose-event", cbExposeCanvas, NULL)


  ##-----------------------
  ##   スケールバーの作成
  ##-----------------------
  ##-----for 'GFI'-----
  hbox.GFI <- gtkHBoxNew (FALSE, 5)
  gtkBoxPackStart (vbox, hbox.GFI, FALSE, FALSE, 0)

  ##GFIの値表示
  alignment.GFI <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.GFI, TRUE, FALSE, 0)

  label.GFI <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.GFI, label.GFI)
  gtkWidgetSetSizeRequest (label.GFI, 130, 50)
  gSignalConnect (label.GFI, "expose-event", cbExposeLabelGFI, NULL)

  #gtkBoxPackStart (hbox.GFI, scale.lambda, TRUE, TRUE, 0)


  ##AGFIの値表示
#  vbox.AGFI <- gtkVBoxNew (FALSE, 5)
#  gtkBoxPackStart (vbox.AGFI, hbox.GFI, FALSE, FALSE, 0)
  alignment.AGFI <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.AGFI, TRUE, FALSE, 0)

  label.AGFI <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.AGFI, label.AGFI)
  gtkWidgetSetSizeRequest (label.AGFI, 130, 50)
  gSignalConnect (label.AGFI, "expose-event", cbExposeLabelAGFI, NULL)
  

  ##AICの値表示
#  vbox.AIC <- gtkVBoxNew (FALSE, 5)
#  gtkBoxPackStart (vbox.AIC, hbox.GFI, FALSE, FALSE, 0)
  alignment.AIC <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.AIC, TRUE, FALSE, 0)

  label.AIC <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.AIC, label.AIC)
  gtkWidgetSetSizeRequest (label.AIC, 130, 50)
  gSignalConnect (label.AIC, "expose-event", cbExposeLabelAIC, NULL)


  ##BICの値表示
#  vbox.BIC <- gtkVBoxNew (FALSE, 5)
#  gtkBoxPackStart (vbox.BIC, hbox.GFI, FALSE, FALSE, 0)
  alignment.BIC <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.BIC, TRUE, FALSE, 0)

  label.BIC <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.BIC, label.BIC)
  gtkWidgetSetSizeRequest (label.BIC, 130, 50)
  gSignalConnect (label.BIC, "expose-event", cbExposeLabelBIC, NULL)


  ##CAICの値表示
#  vbox.CAIC <- gtkVBoxNew (FALSE, 5)
#  gtkBoxPackStart (vbox.CAIC, hbox.GFI, FALSE, FALSE, 0)
  alignment.CAIC <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.CAIC, TRUE, FALSE, 0)

  label.CAIC <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.CAIC, label.CAIC)
  gtkWidgetSetSizeRequest (label.CAIC, 130, 50)
  gSignalConnect (label.CAIC, "expose-event", cbExposeLabelCAIC, NULL)


  ##-----------------------
  ##   スケールバーの作成
  ##-----------------------
  ##-----for 'rho'-----
  hbox.lambda <- gtkHBoxNew (FALSE, 5)
  gtkBoxPackStart (vbox, hbox.lambda, FALSE, FALSE, 0)

  scale.lambda <- gtkHScaleNewWithRange (Min.lambda, Max.lambda, Step.lambda)
  gtkScaleSetDigits (scale.lambda, 2)
  gtkScaleSetDrawValue (scale.lambda, FALSE) ##バーの下に値の表示の有無

  ##rhoの値表示
  alignment.lambda <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.lambda, alignment.lambda, FALSE, FALSE, 0)

  label.lambda <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.lambda, label.lambda)
  gtkWidgetSetSizeRequest (label.lambda, 120, 50)
  #gtkWidgetSetSizeRequest (label.lambda, 100, 50)
  gSignalConnect (label.lambda, "expose-event", cbExposeLabelLambda, NULL)

  gtkBoxPackStart (hbox.lambda, scale.lambda, TRUE, TRUE, 0)



  ##-----for 'gamma'-----
  hbox.gamma <- gtkHBoxNew (FALSE, 5)
  gtkBoxPackStart (vbox, hbox.gamma, FALSE, FALSE, 0)

  scale.gamma <- gtkHScaleNewWithRange (Min.gamma, Max.gamma, Step.gamma)
  gtkScaleSetDigits (scale.gamma, 2)
  gtkScaleSetDrawValue (scale.gamma, FALSE) ##バーの下に値の表示の有無

  ##gammaの値表示
  alignment.gamma <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.gamma, alignment.gamma, FALSE, FALSE, 0)

  label.gamma <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.gamma, label.gamma)
  gtkWidgetSetSizeRequest (label.gamma, 120, 50)
  #gtkWidgetSetSizeRequest (label.gamma, 100, 50)
  gSignalConnect (label.gamma, "expose-event", cbExposeLabelGamma, NULL)

  gtkBoxPackStart (hbox.gamma, scale.gamma, TRUE, TRUE, 0)


  ##value-changedシグナルの発生
  gSignalConnect (scale.lambda, "value_changed", cbValueChangedLambda, canvas)
  gSignalConnect (scale.lambda, "value_changed", cbValueChangedLabelLambda, list(label.lambda, label.GFI, label.AGFI, label.AIC, label.BIC, label.CAIC))
  gSignalConnect (scale.gamma, "value_changed", cbValueChangedGamma, canvas)
  gSignalConnect (scale.gamma, "value_changed", cbValueChangedLabelGamma, list(label.gamma, label.lambda, label.GFI, label.AGFI, label.AIC, label.BIC, label.CAIC))
}


##-----事前にinfo.fancに含めるべき変数-----##

##*Rad.Ellipse: 因子の楕円の長径
##  -> 短径は長径から計算している
##*Len.Rec: 変数の正方形の辺の長さ
##*Window.Height: ウインドウの高さ
##  -> ウインドウの幅は現在は変数及び因子の数から計算している
##*num.lambda: Lambdaの初期番号（変数Lambdaの内の序数）
##*L: Lambdaごとのパラメータ（因子負荷量）を3次元配列で与える．
##  -> 配列の次元は（変数，因子，lambda）に対応している．
##*lambdas, *gammas: lambda, gammasの値をベクトル型で与える.
##*lambda.current, *gamma.current: lambda, gammaの現在値
##*N.var, *N.fac, *N.lambda, *N.gamma: 変数の数，因子数，及びtuning parameterの数
##  -> Lから計算する


##plot.fanc
##gammaの値を変化できるように修正
plot.fanc <- function(x, Window.Height=500, ...){
  if(Window.Height<250 || Window.Height>2000) stop("'Window.Height' must be in [250,2000].")
  if(nchar(system.file(package="RGtk2")) == 0) stop("Package 'RGtk2' is required to plot the solution path.")
  require(RGtk2, quietly=TRUE)
	##if(sum(ls()=="info.fanc")!=0) stop('Object "info.fanc" must be deleted')

	##因子パターン行列の取りだし
  L <- x$loadings
  lambdas <- x$rho
  gammas <- x$gamma
  GFIs <- x$GFI
  AGFIs <- x$AGFI
  AICs <- x$AIC
  BICs <- x$BIC
  CAICs <- x$CAIC
  info.fanc <<- list("Rad.Ellipse"=50, "Len.Rec"=50, "Window.Height"=Window.Height,
                     "N.var"=NULL, "N.fac"=NULL, "N.lambda"=NULL,
                    "L"=NULL, "lambdas"=NULL, "num.lambda"=1, "num.gamma"=1,"num.GFI"=1)
	#if(x$factors==1) L <- array(L,dim=c(nrow(x$rho),1,ncol(x$rho)))
	#assign("info,fanc", info.fanc, envir=infofanc)
  info.fanc$lambda.current <<- lambdas[1,1]
  info.fanc$gamma.current <<- gammas[1]
  info.fanc$GFI.current <<- GFIs[1]
  info.fanc$AGFI.current <<- AGFIs[1]
  info.fanc$AIC.current <<- AICs[1]
  info.fanc$BIC.current <<- BICs[1]
  info.fanc$CAIC.current <<- CAICs[1]
	info.fanc$L <<- L
	info.fanc$lambdas <<- lambdas
	info.fanc$gammas <<- gammas
	info.fanc$GFIs <<- GFIs
	info.fanc$AGFIs <<- AGFIs
	info.fanc$AICs <<- AICs
	info.fanc$BICs <<- BICs
	info.fanc$CAICs <<- CAICs
	info.fanc$N.var <<- dim(info.fanc$L)[1]
	info.fanc$N.fac <<- dim(info.fanc$L)[2]
	info.fanc$N.lambda <<- dim(info.fanc$L)[3]
	info.fanc$N.gamma <<- dim(info.fanc$L)[4]
	info.fanc$Window.Width <<- max(((info.fanc$N.var+1) * info.fanc$Len.Rec + (info.fanc$N.var+2) * info.fanc$Len.Rec / 2), ((info.fanc$N.fac) * 1.5 * info.fanc$Rad.Ellipse),650)


	MakeInterface ("Factor analysis with MC+")
}
