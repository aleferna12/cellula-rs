use crate::pond::Pond;

pub trait Disperser {
    type Error;
    fn disperse(&mut self, source: &mut Pond, target: &mut Pond) -> Result<(), Self::Error>;
}

pub struct SelectiveDispersion {
}

impl SelectiveDispersion {

}

impl Disperser for SelectiveDispersion {
    type Error = ();

    fn disperse(&mut self, source: &mut Pond, target: &mut Pond) -> Result<(), Self::Error> {

    }
}