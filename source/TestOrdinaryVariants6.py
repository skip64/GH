import OrdinaryVariants
import OrdinaryGraphComplex

l = 10
for v in range(20):
    V = OrdinaryVariants.OrdinaryGVSFull(v,l,False)
    W = OrdinaryVariants.OrdinaryGVSTriconnected(v,l,False)
    U = OrdinaryGraphComplex.OrdinaryGVS(v,l,False)
    vv = V.get_dimension()
    ww=W.get_dimension()
    uu=U.get_dimension()
    # print(l,v,":",vv,ww,uu)
    if vv != 0:
        # print("&", v, "&" , vv, "&", uu, "&", "Percentage:",100*ww/vv, 100*uu/vv)
        print(F"& {v} & {vv} & {uu} ({100*uu//vv}\\%) & {ww} ({100*ww//vv}\\%) \\\\")