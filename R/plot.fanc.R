##��������������������������������������������������������
##  Factor analysis with lasso�ɂ�����Gtk�I�u�W�F�N�g
##
##  �t�@�C�����Fplot.fanc.R
##  �t�@�C�����e�F
##  �쐬�ҁFYAMAMOTO, Michio
##  �쐬���F2012�N02��02��
##  �ŏI�X�V���F2012�N03��19��
##  �R�����g�F"X", "F"�̕\�����폜�i031012�j
##           �Ƃ肠�����g����`�ɂ����i031012�j
##           �}�C�i�X�̕��חʂ͐Ԃŕ\�����ƂƂ����i031012�j
##           gamma�̃o�[��ǉ��C���q������̏ꍇ�̃v���b�g�C���i031912�j
##��������������������������������������������������������

##-----�R�[���o�b�N�֐���`-----##
##expose-event�V�O�i���ɑ΂���R�[���o�b�N�֐�
##�}�`�`��֐�
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
  ##   ���q�̑ȉ~��`��
  ##-------------------
  Rem.N.fac <- N.fac %% 2
  cairoTranslate (cr, Window.Width / 2, Rad.Ellipse + 20) ##���_�̈ړ�
  cairoTranslate (cr.t, Window.Width / 2, Rad.Ellipse + 20) ##���_�̈ړ�

  cairoSetLineWidth (cr, 2.5)

  cairoScale (cr, 1.0, 0.5) ##�ȉ~�̃X�P�[��

  if (Rem.N.fac == 0) { ##���q���������̏ꍇ
    for (i in 1:HN.fac) {
      ii <- i - 1

      ##�ȉ~�̕`��
      cairoArc (cr, -Rad.Ellipse * 5 / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 0,
                Rad.Ellipse, 0.0, 2.0 * pi)
      cairoStroke (cr)
      cairoArc (cr,  Rad.Ellipse * 5 / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 0,
                Rad.Ellipse, 0.0, 2.0 * pi)
      cairoStroke (cr)

      cairoSetFontSize (cr.t, 20)

      ##���x���̕`��
      text <- sprintf ("f%d", HN.fac - ii)
      cairoMoveTo (cr.t, -Rad.Ellipse * 5 / 4 - Rad.Ellipse / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 5)
      cairoShowText (cr.t, text)

      text <- sprintf ("f%d", HN.fac + 1 + ii)
      cairoMoveTo (cr.t,  Rad.Ellipse * 5 / 4 - Rad.Ellipse / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 5)
      cairoShowText (cr.t, text)
    }
  }
  else if (Rem.N.fac != 0) { ##���q������̏ꍇ
    cairoArc (cr, 0.0, 0.0, Rad.Ellipse, 0.0, 2.0 * pi)
    cairoStroke (cr)
    cairoSetFontSize (cr.t, 20)
    text <- sprintf ("f%d", floor(HN.fac) + 1)
    cairoMoveTo (cr.t, 0 - Rad.Ellipse / 4, 5)
    cairoShowText (cr.t, text)

    if (floor(HN.fac) != 0) {
      for (i in 1:floor(HN.fac)) {
        ii <- i - 1

        ##�ȉ~�̕`��
        cairoArc (cr, -Rad.Ellipse * 10 / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 0,
                  Rad.Ellipse, 0.0, 2.0 * pi)
        cairoStroke (cr)
        cairoArc (cr,  Rad.Ellipse * 10 / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 0,
                  Rad.Ellipse, 0.0, 2.0 * pi)
        cairoStroke (cr)

        cairoSetFontSize (cr.t, 20)

        ##���x���̕`��
        text <- sprintf ("f%d", floor(HN.fac) - ii)
        cairoMoveTo (cr.t, -Rad.Ellipse * 10 / 4 - Rad.Ellipse / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 5)
        cairoShowText (cr.t, text)

        text <- sprintf ("f%d", floor(HN.fac) + ii + 2)
        cairoMoveTo (cr.t,  Rad.Ellipse * 10 / 4 - Rad.Ellipse / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 5)
        cairoShowText (cr.t, text)
      }
    }
  } ##���q�̑ȉ~�̕`��I��

  ##����܂ł̍��W�n�̕ϊ������ɖ߂�
  cairoGetMatrix (cr, Mat)
  cairoMatrixInvert (Mat)
  cairoTransform (cr, Mat)
  cairoGetMatrix (cr.t, Mat)
  cairoMatrixInvert (Mat)
  cairoTransform (cr.t, Mat)


  ##-------------------
  ##  �ϐ��̎l�p��`��
  ##-------------------
  ##�ϐ��̐��̋��ŏꍇ�킯
  Rem.N.var <- N.var %% 2
  cairoTranslate (  cr, Window.Width / 2, Window.Height - Len.Rec * 5 / 2)
  cairoTranslate (cr.t, Window.Width / 2, Window.Height - Len.Rec * 5 / 2)

  cairoSetLineWidth (cr, 2.5)

  if (Rem.N.var == 0) {##�ϐ��̐��������̏ꍇ
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
  else if (Rem.N.var != 0) { ##�ϐ��̐�����̏ꍇ
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

  ##����܂ł̍��W�n�̕ϊ������ɖ߂�
  cairoGetMatrix (cr, Mat)
  cairoMatrixInvert (Mat)
  cairoTransform (cr, Mat)
  cairoGetMatrix (cr.t, Mat)
  cairoMatrixInvert (Mat)
  cairoTransform (cr.t, Mat)



  ##---------------------------
  ##  ���q�ƕϐ������Ԑ���`��
  ##---------------------------
  if (Rem.N.fac == 0) { ##���q���������̏ꍇ
    cr.Fac.X <- -Rad.Ellipse * 5 / 4 - (HN.fac - 1) * 2 * Rad.Ellipse * 5 / 4 + Window.Width / 2
    cr.Fac.Y <- Rad.Ellipse + 20 + Rad.Ellipse / 2
  }
  else if (Rem.N.fac != 0) { ##���q������̏ꍇ
    cr.Fac.X <- -Rad.Ellipse * 10 / 4 - (floor(HN.fac) - 1) * 2 * Rad.Ellipse * 5 / 4 + Window.Width / 2
    cr.Fac.Y <- Rad.Ellipse + 20 + Rad.Ellipse / 2
  }

  if (Rem.N.var == 0) { ##�ϐ��̐��������̏ꍇ
    cr.Var.X <- -Len.Rec * 5 / 4 - (HN.var - 1) * Len.Rec * 3 / 2 + Window.Width / 2 + Len.Rec / 2
    cr.Var.Y <- Window.Height - Len.Rec * 5 / 2
  }
  else if (Rem.N.var != 0) { ##�ϐ��̐�����̏ꍇ
    cr.Var.X <- -Len.Rec * 2 - (floor(HN.var) - 1) * Len.Rec * 3 / 2 + Window.Width / 2 + Len.Rec / 2
    cr.Var.Y <- Window.Height - Len.Rec * 5 / 2
  }

  ##�S�Ẵp�X���X�ɕ`�悷��
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



##GtkWidget label��expose-event�V�O�i���ɑ΂���R�[���o�b�N�֐�
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


##GtkScale�E�B�W�F�b�g�i"canvas"�j��value-changed�V�O�i���ɑ΂���R�[���o�b�N�֐�
##-----for 'rho'-----
cbValueChangedLambda <- function (gtk.scale, data)
{
  if (gtkRangeGetValue (gtk.scale) == 0) {
    info.fanc$num.lambda <<- 1
  }
  else {
    info.fanc$num.lambda <<- ceiling(gtkRangeGetValue (gtk.scale) * info.fanc$N.lambda) ##�덷���C��
  }
  gtkWidgetQueueDraw (data) ##canvas�ɑ΂���expose-event�V�O�i�����o��
}

##-----for 'gamma'-----
cbValueChangedGamma <- function (gtk.scale, data)
{
  if (gtkRangeGetValue (gtk.scale) == 0) {
    info.fanc$num.gamma <<- 1
  }
  else {
    info.fanc$num.gamma <<- ceiling(gtkRangeGetValue (gtk.scale) * info.fanc$N.gamma) ##�덷���C��
  }

  info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]

  gtkWidgetQueueDraw (data) ##canvas�ɑ΂���expose-event�V�O�i�����o��
}



##GtkScale�E�B�W�F�b�g�i"label"�j��value-changed�V�O�i���ɑ΂���R�[���o�b�N�֐�
##-----for 'rho'-----
cbValueChangedLabelLambda <- function (gtk.scale, data)
{
  #if (gtkRangeGetValue (gtk.scale) == 0) {
  #  info.fanc$lambda.current <<- info.fanc$lambdas[1, info.fanc$num.gamma]
  #}
  #else {
  #  info.fanc$lambda.current <<- info.fanc$lambdas[ceiling(gtkRangeGetValue (gtk.scale) * info.fanc$N.lambda), info.fanc$num.gamma] ##�덷���C��
  #}
  info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$GFI.current <<- info.fanc$GFIs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$AGFI.current <<- info.fanc$AGFIs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$AIC.current <<- info.fanc$AICs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$BIC.current <<- info.fanc$BICs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$CAIC.current <<- info.fanc$CAICs[info.fanc$num.lambda, info.fanc$num.gamma]

  gtkWidgetQueueDraw (data[[1]]) ##label.lambda�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[2]]) ##label.GFI�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[3]]) ##label.AGFI�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[4]]) ##label.AIC�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[5]]) ##label.BIC�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[6]]) ##label.CAIC�ɑ΂���expose-event�V�O�i�����o��
}

##-----for 'gamma'-----
cbValueChangedLabelGamma <- function (gtk.scale, data)
{
  if (gtkRangeGetValue (gtk.scale) == 0) {
    info.fanc$gamma.current <<- info.fanc$gammas[1]
  }
  else {
    info.fanc$gamma.current <<- info.fanc$gammas[ceiling(gtkRangeGetValue (gtk.scale) * info.fanc$N.gamma)] ##�덷���C��
  }

  info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$GFI.current <<- info.fanc$GFIs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$AGFI.current <<- info.fanc$AGFIs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$AIC.current <<- info.fanc$AICs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$BIC.current <<- info.fanc$BICs[info.fanc$num.lambda, info.fanc$num.gamma]
  info.fanc$CAIC.current <<- info.fanc$CAICs[info.fanc$num.lambda, info.fanc$num.gamma]

  gtkWidgetQueueDraw (data[[1]]) ##label.gamma�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[2]]) ##label.lambda�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[3]]) ##label.GFI�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[4]]) ##label.AGFI�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[5]]) ##label.AIC�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[6]]) ##label.BIC�ɑ΂���expose-event�V�O�i�����o��
  gtkWidgetQueueDraw (data[[7]]) ##label.CAIC�ɑ΂���expose-event�V�O�i�����o��
}



##-----�`�悷�邽�߂̃E�B���h�E�{��-----##
MakeInterface <- function (gchar.title)
{
  ##"info" includes 'L' and 'lambdas'.
  ##'L' is a p*m*r array of pattern matrices
  ##'lambdas' is a r*1 vector of tuning parameter

  lambdas <- info.fanc$lambdas
  gammas <- info.fanc$gammas

  ##GTK+��ł�lambda, gamma�͈̔͂�(0, 1)�Ɍ��肷��
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
##  gSignalConnect (window, "destroy", gtkMainQuit, NULL) ##�Ȃ��Ă��ǂ���������Ȃ�


  ##-------------------
  ##    �p�X�}�̕`��
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
  ##   �X�P�[���o�[�̍쐬
  ##-----------------------
  ##-----for 'GFI'-----
  hbox.GFI <- gtkHBoxNew (FALSE, 5)
  gtkBoxPackStart (vbox, hbox.GFI, FALSE, FALSE, 0)

  ##GFI�̒l�\��
  alignment.GFI <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.GFI, TRUE, FALSE, 0)

  label.GFI <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.GFI, label.GFI)
  gtkWidgetSetSizeRequest (label.GFI, 130, 50)
  gSignalConnect (label.GFI, "expose-event", cbExposeLabelGFI, NULL)

  #gtkBoxPackStart (hbox.GFI, scale.lambda, TRUE, TRUE, 0)


  ##AGFI�̒l�\��
#  vbox.AGFI <- gtkVBoxNew (FALSE, 5)
#  gtkBoxPackStart (vbox.AGFI, hbox.GFI, FALSE, FALSE, 0)
  alignment.AGFI <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.AGFI, TRUE, FALSE, 0)

  label.AGFI <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.AGFI, label.AGFI)
  gtkWidgetSetSizeRequest (label.AGFI, 130, 50)
  gSignalConnect (label.AGFI, "expose-event", cbExposeLabelAGFI, NULL)
  

  ##AIC�̒l�\��
#  vbox.AIC <- gtkVBoxNew (FALSE, 5)
#  gtkBoxPackStart (vbox.AIC, hbox.GFI, FALSE, FALSE, 0)
  alignment.AIC <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.AIC, TRUE, FALSE, 0)

  label.AIC <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.AIC, label.AIC)
  gtkWidgetSetSizeRequest (label.AIC, 130, 50)
  gSignalConnect (label.AIC, "expose-event", cbExposeLabelAIC, NULL)


  ##BIC�̒l�\��
#  vbox.BIC <- gtkVBoxNew (FALSE, 5)
#  gtkBoxPackStart (vbox.BIC, hbox.GFI, FALSE, FALSE, 0)
  alignment.BIC <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.BIC, TRUE, FALSE, 0)

  label.BIC <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.BIC, label.BIC)
  gtkWidgetSetSizeRequest (label.BIC, 130, 50)
  gSignalConnect (label.BIC, "expose-event", cbExposeLabelBIC, NULL)


  ##CAIC�̒l�\��
#  vbox.CAIC <- gtkVBoxNew (FALSE, 5)
#  gtkBoxPackStart (vbox.CAIC, hbox.GFI, FALSE, FALSE, 0)
  alignment.CAIC <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.GFI, alignment.CAIC, TRUE, FALSE, 0)

  label.CAIC <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.CAIC, label.CAIC)
  gtkWidgetSetSizeRequest (label.CAIC, 130, 50)
  gSignalConnect (label.CAIC, "expose-event", cbExposeLabelCAIC, NULL)


  ##-----------------------
  ##   �X�P�[���o�[�̍쐬
  ##-----------------------
  ##-----for 'rho'-----
  hbox.lambda <- gtkHBoxNew (FALSE, 5)
  gtkBoxPackStart (vbox, hbox.lambda, FALSE, FALSE, 0)

  scale.lambda <- gtkHScaleNewWithRange (Min.lambda, Max.lambda, Step.lambda)
  gtkScaleSetDigits (scale.lambda, 2)
  gtkScaleSetDrawValue (scale.lambda, FALSE) ##�o�[�̉��ɒl�̕\���̗L��

  ##rho�̒l�\��
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
  gtkScaleSetDrawValue (scale.gamma, FALSE) ##�o�[�̉��ɒl�̕\���̗L��

  ##gamma�̒l�\��
  alignment.gamma <- gtkAlignmentNew (0.5, 0.5, 0, 0)
  gtkBoxPackStart (hbox.gamma, alignment.gamma, FALSE, FALSE, 0)

  label.gamma <- gtkDrawingAreaNew ()
  gtkContainerAdd (alignment.gamma, label.gamma)
  gtkWidgetSetSizeRequest (label.gamma, 120, 50)
  #gtkWidgetSetSizeRequest (label.gamma, 100, 50)
  gSignalConnect (label.gamma, "expose-event", cbExposeLabelGamma, NULL)

  gtkBoxPackStart (hbox.gamma, scale.gamma, TRUE, TRUE, 0)


  ##value-changed�V�O�i���̔���
  gSignalConnect (scale.lambda, "value_changed", cbValueChangedLambda, canvas)
  gSignalConnect (scale.lambda, "value_changed", cbValueChangedLabelLambda, list(label.lambda, label.GFI, label.AGFI, label.AIC, label.BIC, label.CAIC))
  gSignalConnect (scale.gamma, "value_changed", cbValueChangedGamma, canvas)
  gSignalConnect (scale.gamma, "value_changed", cbValueChangedLabelGamma, list(label.gamma, label.lambda, label.GFI, label.AGFI, label.AIC, label.BIC, label.CAIC))
}


##-----���O��info.fanc�Ɋ܂߂�ׂ��ϐ�-----##

##*Rad.Ellipse: ���q�̑ȉ~�̒��a
##  -> �Z�a�͒��a����v�Z���Ă���
##*Len.Rec: �ϐ��̐����`�̕ӂ̒���
##*Window.Height: �E�C���h�E�̍���
##  -> �E�C���h�E�̕��͌��݂͕ϐ��y�ш��q�̐�����v�Z���Ă���
##*num.lambda: Lambda�̏����ԍ��i�ϐ�Lambda�̓��̏����j
##*L: Lambda���Ƃ̃p�����[�^�i���q���חʁj��3�����z��ŗ^����D
##  -> �z��̎����́i�ϐ��C���q�Clambda�j�ɑΉ����Ă���D
##*lambdas, *gammas: lambda, gammas�̒l���x�N�g���^�ŗ^����.
##*lambda.current, *gamma.current: lambda, gamma�̌��ݒl
##*N.var, *N.fac, *N.lambda, *N.gamma: �ϐ��̐��C���q���C�y��tuning parameter�̐�
##  -> L����v�Z����


##plot.fanc
##gamma�̒l��ω��ł���悤�ɏC��
plot.fanc <- function(x, Window.Height=500, ...){
  if(Window.Height<250 || Window.Height>2000) stop("'Window.Height' must be in [250,2000].")
  if(nchar(system.file(package="RGtk2")) == 0) stop("Package 'RGtk2' is required to plot the solution path.")
  require(RGtk2, quietly=TRUE)
	##if(sum(ls()=="info.fanc")!=0) stop('Object "info.fanc" must be deleted')

	##���q�p�^�[���s��̎�肾��
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