#[repr(C)]
pub struct KinectHandle {
    _data: (),
    _marker:
        core::marker::PhantomData<(*mut u8, core::marker::PhantomPinned)>,
}

unsafe extern "C" {
    pub fn kinect_create(rgb: bool, depth: bool) -> *mut KinectHandle;
    pub fn kinect_listen_frame(h: *mut KinectHandle, ms: i32) -> bool;
    pub fn kinect_release_frame(h: *mut KinectHandle);
    pub fn kinect_destroy(h: *const KinectHandle);
    pub fn kinect_depth(h: *mut KinectHandle) -> *const f32;
}