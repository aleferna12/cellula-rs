use bitvec::slice::BitSlice;

pub struct BitAdhesion;

impl BitAdhesion {
    pub fn contact_energy(lig: &BitSlice<u64>, rec: &BitSlice<u64>) -> u16 {
        // There is probably a faster way but this is only called during initialization anyway
        lig.iter().by_vals().zip(rec.iter().by_vals()).filter(|(l, r)| l ^ r).count() as u16
    }
}