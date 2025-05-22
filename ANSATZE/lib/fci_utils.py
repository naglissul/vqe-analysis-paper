from pyscf import gto, scf, ao2mo, fci

def get_fci_h2(distance: float = 0.735) -> float:
    mol = gto.Mole()
    mol.atom = f'H 0 0 0; H 0 0 {distance}'
    mol.basis = 'sto-3g'
    mol.spin = 0
    mol.charge = 0
    mol.build()

    mf = scf.RHF(mol)
    mf.kernel()

    norb = mf.mo_coeff.shape[1]
    nelec = mol.nelec

    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.kernel(mol, mf.mo_coeff)

    cisolver = fci.FCI(mol, mf.mo_coeff)
    E_fci, _ = cisolver.kernel(h1, eri, norb, nelec)

    return E_fci

def get_fci_lih(distance: float = 1.6) -> float:
    mol = gto.Mole()
    mol.atom=f'Li 0 0 0; H 0 0 {distance}'
    mol.basis = 'sto-3g'
    mol.spin = 0
    mol.charge = 0
    mol.build()

    mf = scf.RHF(mol)
    mf.kernel()

    norb = mf.mo_coeff.shape[1]
    nelec = mol.nelec

    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.kernel(mol, mf.mo_coeff)

    cisolver = fci.FCI(mol, mf.mo_coeff)
    E_fci, _ = cisolver.kernel(h1, eri, norb, nelec)

    return E_fci