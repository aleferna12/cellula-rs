#[repr(C)]
pub struct KinectHandle {
    _data: (),
    _marker:
        core::marker::PhantomData<(*mut u8, core::marker::PhantomPinned)>,
}

unsafe extern "C" {
    pub fn kinect_create() -> *mut KinectHandle;
    pub fn kinect_next_depth(h: *mut KinectHandle, ms: i32) -> *const f32;
    pub fn kinect_release_frame(h: *mut KinectHandle);
    pub fn kinect_destroy(h: *const KinectHandle);
}