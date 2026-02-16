import numpy as np

class SFS_1D_XXZ:
    def __init__(self, nmo, delta):
        """
        OBC, p runs 1..N-1 in 1-based language.
        0-based ↔ 1-based:
            site i (0-based) ↔ site (i+1) (1-based)
        Ghost convention: N_0 = 2.
        """
        self.nmo = nmo
        self.delta = delta

        # ---------- h100, h001 ----------
        # From (1/4) * N_{p-1} P_p^(†) at p=1 -> (1/4)*N0 * P_1^(†) with N0=2 => 1/2
        self.h100 = np.zeros(nmo, dtype=float)
        self.h001 = np.zeros(nmo, dtype=float)
        if nmo >= 1:
            self.h100[0] = 0.50
            self.h001[0] = 0.50

        # ---------- h010 ----------
        # Build the Δ-dependent one-body N terms correctly for all nmo,
        # including the special small-N cases.
        self.h010 = np.zeros(nmo, dtype=float)

        if nmo == 1:
            # No p in 1..N-1, so no Δ one-body or constant from those sums.
            pass

        elif nmo == 2:
            # N=2 => p=1 only.
            # One-body:  -Δ/4 N_{p+1} gives -Δ/4 N2 (index 1)
            #            +Δ/4 N_{p-1}N_{p+1} gives +(Δ/4) N0 N2 = +Δ/2 N2
            # => h010[1] = +Δ/4
            self.h010[0] = 0.0
            self.h010[1] = +0.25 * delta

        elif nmo == 3:
            # N=3 => p=1,2.
            # Result: h010 = [-Δ/4, +Δ/4, -Δ/4]
            self.h010[0] = -0.25 * delta
            self.h010[1] = +0.25 * delta
            self.h010[2] = -0.25 * delta

        else:
            # nmo >= 4 (the pattern we discussed for longer chains):
            # index 0 : -Δ/4
            # index 1 : 0 (cancellation)
            # index 2..nmo-3 : -Δ/2
            # index nmo-2, nmo-1 : -Δ/4
            self.h010[0] = -0.25 * delta
            self.h010[1] = 0.0
            for i in range(2, nmo - 2):
                self.h010[i] = -0.5 * delta
            self.h010[nmo - 2] = -0.25 * delta
            self.h010[nmo - 1] = -0.25 * delta

        # ---------- h110 ----------
        self.h110 = np.zeros((nmo, nmo), dtype=float)
        if nmo >= 2:
            self.h110[0, 1] = -0.25
        for p in range(1, nmo - 1):     # p = 1 .. nmo-2
            self.h110[p, p - 1] += 0.25
            self.h110[p, p + 1] += 0.25

        # ---------- h011 ----------
        self.h011 = np.zeros((nmo, nmo), dtype=float)
        if nmo >= 2:
            self.h011[1, 0] = -0.25
        for p in range(1, nmo - 1):
            self.h011[p - 1, p] += 0.25
            self.h011[p + 1, p] += 0.25

        # ---------- h120 ----------
        self.h120 = np.zeros((nmo, nmo, nmo), dtype=float)
        for p in range(1, nmo - 1):
            q = p - 1
            r = p + 1
            self.h120[p, q, r] += -0.125
            self.h120[p, r, q] += -0.125

        # ---------- h021 ----------
        self.h021 = np.zeros((nmo, nmo, nmo), dtype=float)
        for p in range(1, nmo - 1):
            left = p - 1
            right = p + 1
            self.h021[left, right, p] += -0.125
            self.h021[right, left, p] += -0.125

        # ---------- h020 ----------
        self.h020 = np.zeros((nmo, nmo), dtype=float)
        for p in range(1, nmo - 1):
            i = p - 1
            j = p + 1
            self.h020[i, j] += 0.125 * delta
            self.h020[j, i] += 0.125 * delta

        # ---------- h000 ----------
        # Your formula works for nmo >= 2 (gives 0 at nmo=3, -Δ/4 at nmo=2).
        # For nmo=1 it should be 0 because there is no p-sum at all.
        if nmo == 1:
            self.h000 = 0.0
        else:
            self.h000 = 0.25 * delta * (nmo - 3)
