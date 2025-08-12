use crate::ecology::disperser::DispersionEvent;

pub trait Transporter {
    type Error;
    fn transport(&mut self, event: DispersionEvent) -> Result<(), Self::Error>;
}