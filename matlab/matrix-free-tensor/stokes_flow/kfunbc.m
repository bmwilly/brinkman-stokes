function w = kfunbc(u, xy, xyp, mv, bound, ae, bxe, bye)

wa = afunbc(u, xy, xyp, mv, bound, ae);
wb = bfunbc(u, xy, xyp, mv, bound, bxe, bye);
wbt = btfunbc(u, xy, xyp, mv, bound, bxe, bye);
w = wa + wb + wbt;

end