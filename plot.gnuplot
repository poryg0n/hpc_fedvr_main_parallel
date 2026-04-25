set term pngcairo

set output "psi.png"
plot "final_psi.dat"    using 2:3 w l title "Re({/Symbol Y})(x,t)", \
     "final_psi_ex.dat" using 2:3 w l title "Re({/Symbol Y})_{ex}(x,t)"


set output "phi.png"
plot "final_phi.dat"    using 2:3 w l title "Re({/Symbol F})(x,t)", \
     "final_phi_ex.dat" using 2:3 w l title "Re({/Symbol F})_{ex}(x,t)", \
